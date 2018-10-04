"""
OFFLINE PHASE aerostruct
This program aims to use the POD-ROM to create a reduced basis able to be used in a real-time application.
It is applied to an aeroelastic problem. A Greedy algorithm has been used to introduced a POD based on 
the snapshots method combined with the Singular Value Decomposition (SVD).
Using version 3.0.6 of GMSH.

July 2018
@author: ochandre
"""

import numpy as np
from scipy import linalg
from pyDOE import lhs
import functions_Offline as foff
import timeit
from sklearn.externals import joblib
from gaussian_process import GaussianProcess


##PARAMETERS BOUNDS: h_skins, h_ribs, h_spars_le, h_spars_te, b, S
# h_* = thickness of different parts (m)
# b = wing span (m)
# S = wing surface (m^2)
pMin_i = np.array([0.028,0.009,0.017,0.009,64.75,443.56])
pMax_i = np.array([0.031,0.015,0.025,0.015,75.65,500.21])
pMin = np.zeros((6,1))
pMin[:,0]= pMin_i[:]
pMax = np.zeros((6,1))
pMax[:,0]= pMax_i[:]
np.save("../results/Offline/pMin",pMin)
np.save("../results/Offline/pMax",pMax)

##OTHERS PARAMETERS OF INTEREST
# alpha = angle of attack (°)
# M_inf = Mach number
# E = Young modulus (Pa)
# nu = Poisson's ratio
# BF = distance of the fuselage over the wing (m)
# phi_d = sweep angle (°)
# diedre_d = dihedral angle (°)
# n_ribs_1 = number of ribs (first sector)
# n_ribs_2 = number of ribs (second sector)
alpha = 2.5 
M_inf = 0.85
E = 70.*1e9
nu = 0.3
BF = 5.96
phi_d = 30.
diedre_d = 1.5
n_ribs_1 =  8
n_ribs_2 =  12

##FIXED PARAMETERS (Flight Altitude = 11000 meters)
# rho = density at 11000 meters (kg/m^3)
# Pressure = pressure at 11000 meters (Pa)
# adiabatic_i = adiabatic index of air
# a_speed = speed of sound at 11000 meters (m/s)
# phi = sweep angle (rad)
# diedre = dihedral angle (rad)
rho = 0.3629
Pressure = 22552.
adiabiatic_i = 1.4
a_speed = np.sqrt(adiabiatic_i*Pressure/rho)
v_inf = M_inf*a_speed
phi = phi_d*np.pi/180.
diedre = diedre_d*np.pi/180.

##CREATION OF MESH FILES VLM AND FEM
foff.create_mesh_files(n_ribs_1,n_ribs_2)

param_data = np.zeros(9)
param_data[0] = alpha
param_data[1] = v_inf
param_data[2] = rho
param_data[3] = E
param_data[4] = nu
np.save("../results/Offline/param_data",param_data)

##GREEDY ALGORITHM
tic = timeit.default_timer()
# 0) Initialization
nSamples = 3
nIterations = nSamples-1
nCandidates = 5
nParam = len(pMin)
pCandidate = np.dot(pMin,np.ones((1,nCandidates))) + lhs(nParam,nCandidates).T*np.dot((pMax-pMin),np.ones((1,nCandidates)))
indexCand = np.zeros((nCandidates,), dtype=int)
for i in range(nCandidates):
    indexCand[i] = i
pCandMax = np.zeros((nParam,nSamples))
np.save("../results/Offline/pCandidate",pCandidate)
np.save("../results/Offline/nCandidates",nCandidates)
np.save("../results/Offline/nSamples",nSamples)
# 1) Ramdomly select a first sample (in the middle of the domain)
pMiddle = (pMin+pMax)/2
param_data[5] = pMiddle[0,0]
param_data[6] = pMiddle[1,0]
param_data[7] = pMiddle[2,0]
param_data[8] = pMiddle[3,0]
b = pMiddle[4,0]
S = pMiddle[5,0]
pCandMax[:,0] = pMiddle[:,0]
# Calculating the number of nodes and initialization of samples matrices
n_VLM_nodes, n_FEM_nodes, n_gamma_nodes = foff.calcule_nodes(param_data)
uSamples = np.zeros((n_FEM_nodes*6,nSamples)) # POD displacement
gSamples = np.zeros((n_gamma_nodes,nSamples)) # POD circulation
# 2) Solving HDM-based problem
print "Greedy Iteration num.: 0\n",
uMiddle,gMiddle = foff.run_aerostruct(param_data,b,S,phi,diedre,BF,Mach=M_inf)
uSamples[:,0] = uMiddle
gSamples[:,0] = gMiddle
# 3) Building the ROB
uV,_,_ = linalg.svd(uSamples,full_matrices=False)
gV,_,_ = linalg.svd(gSamples,full_matrices=False)
uV_tr = np.transpose(uV)
gV_tr = np.transpose(gV)
# 4) Iteration loop
for iIteration in range(nIterations):
    print "Greedy Iteration num.: ",iIteration+1
    # 5) Solve argmax(residual)
    maxError = 0.
    candidateMaxError = 0
    for iCandidate in indexCand:
        # 5) Loading candidates information
        param_data[5] = pCandidate[0,iCandidate]
        param_data[6] = pCandidate[1,iCandidate]
        param_data[7] = pCandidate[2,iCandidate]
        param_data[8] = pCandidate[3,iCandidate]
        b = pCandidate[4,iCandidate]
        S = pCandidate[5,iCandidate]
        A, B = foff.VLM_params(param_data,b,S,phi,diedre,BF,Mach=M_inf)
        Ar = np.dot(np.dot(gV_tr,A),gV)
        Br = np.dot(gV_tr,B)
        Ar_lu = linalg.lu_factor(Ar)
        gr = linalg.lu_solve(Ar_lu,-Br)
        g_exp = np.dot(gV,gr)
        F_s = foff.get_Fs(g_exp,param_data)
        K, Fs = foff.get_FEM(param_data,F_s)
        Kr = np.dot(np.dot(uV_tr,K),uV)
        Fr = np.dot(uV_tr,Fs)
        q = linalg.solve(Kr,Fr)
        # 5b) Compute error indicator
        Y1 = np.dot(K,np.dot(uV,q))- Fs
        errorIndicator = np.dot(np.transpose(Y1),Y1)
        if (iCandidate == 0 or maxError < errorIndicator):
            candidateMaxError = iCandidate
            maxError = errorIndicator
    # 6) Solving HDM-based problem for the candidate with the highest error indicator
    param_data[5] = pCandidate[0,candidateMaxError]
    param_data[6] = pCandidate[1,candidateMaxError]
    param_data[7] = pCandidate[2,candidateMaxError]
    param_data[8] = pCandidate[3,candidateMaxError]
    b = pCandidate[4,candidateMaxError]
    S = pCandidate[5,candidateMaxError]
    pCandMax[:,iIteration+1] = pCandidate[:,candidateMaxError]
    uIter,gIter = foff.run_aerostruct(param_data,b,S,phi,diedre,BF,Mach=M_inf)
    uSamples[:,iIteration+1] = uIter
    gSamples[:,iIteration+1] = gIter
    # 7) Building the ROB
    uV,_,_ = linalg.svd(uSamples[:,0:iIteration+2],full_matrices=False)
    gV,_,_ = linalg.svd(gSamples[:,0:iIteration+2],full_matrices=False)
    uV_tr = np.transpose(uV)
    gV_tr = np.transpose(gV)
    ##SAVE OFFLINE DATA
    u_name_V = "../results/Offline/uV"+str(iIteration+1)
    g_name_V = "../results/Offline/gV"+str(iIteration+1)
    np.save(u_name_V,uV)
    np.save(g_name_V,gV)
    np.save("../results/Offline/last_iIteration",iIteration)
    ##UPDATING indexCand: The current index must be removed!
    xC = 0
    while True:
        if indexCand[xC] == candidateMaxError:
            break
        xC += 1
    indexCand = np.delete(indexCand,xC)

toc = timeit.default_timer()
print("TOTAL COMPUTATION TIME: "+str((toc-tic)/60.)+" min")

##KRIGING
# Initialization
tic = timeit.default_timer()
U_comp = uSamples
G_comp = gSamples
np.save("../results/Offline/U_comp",U_comp)
np.save("../results/Offline/G_comp",G_comp)
# Weight calculations
a_kri = np.zeros((nSamples,nSamples))
b_kri = np.zeros((nSamples,nSamples))
for i in range(nSamples):
    for j in range(nSamples):
        a_kri[j,i] = np.dot(U_comp[:,j].T,uV[:,i])
        b_kri[j,i] = np.dot(G_comp[:,j].T,gV[:,i])
np.save("../results/Offline/a_kri",a_kri)
np.save("../results/Offline/b_kri",b_kri)
# Trained solution
mean = "constant"
covariance = "squared_exponential"
theta_U = np.array([100000.0]*6)
theta_L = np.array([0.001]*6)
theta_0 = np.array([1.0]*6)
for i in range(nSamples):
    GP_u = GaussianProcess(regr = mean, corr = covariance,theta0 = theta_0,thetaL = theta_L, thetaU = theta_U)
    GP_u.fit(pCandMax.T, a_kri[:,i])
    GP_g = GaussianProcess(regr = mean, corr = covariance,theta0 = theta_0,thetaL = theta_L, thetaU = theta_U)
    GP_g.fit(pCandMax.T, b_kri[:,i])
    joblib.dump(GP_u, "../results/Offline/GP_alpha_"+str(i)+".pkl")
    joblib.dump(GP_g, "../results/Offline/GP_beta_"+str(i)+".pkl")
    
toc = timeit.default_timer()
print("KRIGING COMPUTATION TIME: "+str(toc-tic)+" s")