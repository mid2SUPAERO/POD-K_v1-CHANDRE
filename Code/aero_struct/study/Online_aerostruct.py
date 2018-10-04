"""
ONLINE PHASE aerostruct
This program aims to use the POD-ROM to create a reduced basis able to be used in a real-time application.
It is applied to an aeroelastic problem.
Using version 3.0.6 of GMSH.

July 2018
@author: ochandre
"""

import numpy as np
from sklearn.externals import joblib
import functions_Offline as foff
from mesh_deformation_airbus import mesh_deformation_airbus
from VLM import VLM_study
from FEM import FEM_study
import timeit
import subprocess

##LOADING... (FROM OFFLINE PROCESS)
nSamples = np.load("../results/Offline/nSamples.npy")
param_data = np.load("../results/Offline/param_data.npy")
iIteration = np.load("../results/Offline/last_iIteration.npy")
u_name = "../results/Offline/uV"+str(iIteration+1)+".npy"
g_name = "../results/Offline/gV"+str(iIteration+1)+".npy"
uV = np.load(u_name)
gV = np.load(g_name)

##ONLINE TREATMENT
#Parameters initialisation
tic = timeit.default_timer()
h_skins = 0.028
h_ribs = 0.009
h_spars_le = 0.017
h_spars_te = 0.009
b = 64.75
S = 443.56

x_conf = np.array([h_skins,h_ribs,h_spars_le,h_spars_te,b,S])
pC = np.zeros((6))
pC[0] = param_data[5] = h_skins
pC[1] = param_data[6] = h_ribs
pC[2] = param_data[7] = h_spars_le
pC[3] = param_data[8] = h_spars_te
pC[4] = b
pC[5] = S
# 1) Kriging prediction of the interested point
u_mp = np.zeros((1,nSamples))
u_vp = np.zeros((1,nSamples))
gamma_mp = np.zeros((1,nSamples))
gamma_vp = np.zeros((1,nSamples))
for i in range(nSamples):
    # *Loading Kriging data
    name_a = joblib.load("../results/Offline/GP_alpha_"+str(i)+".pkl")
    name_b = joblib.load("../results/Offline/GP_beta_"+str(i)+".pkl")
    # **Prediction
    u_mp[0,i], u_vp[0,i] = name_a.predict(x_conf.reshape((1,6)), eval_MSE=True)
    gamma_mp[0,i], gamma_vp[0,i] = name_b.predict(x_conf.reshape((1,6)), eval_MSE=True)
u_mean = 0
u_var2 = 0
gamma_mean = 0
gamma_var2 = 0
for i in range(nSamples):
    u_mean += np.dot(u_mp[0,i],uV[:,i])
    u_var2 += np.dot(u_vp[0,i]**2,uV[:,i]**2)
    gamma_mean += np.dot(gamma_mp[0,i],gV[:,i])
    gamma_var2 += np.dot(gamma_vp[0,i]**2,gV[:,i]**2)
# 2) Calculation of u_est and g_est
u_est = u_mean
u_est1 = u_mean+3.0*np.sqrt(u_var2)
u_est2 = u_mean-3.0*np.sqrt(u_var2)
g_est = gamma_mean
g_est1 = gamma_mean+3.0*np.sqrt(gamma_var2)
g_est2 = gamma_mean-3.0*np.sqrt(gamma_var2)
toc = timeit.default_timer()
print("ONLINE COMPUTATION TIME: "+str(toc-tic)+" s")
# 3) Publication of the results: u_est, g_est
a_new = mesh_deformation_airbus('../mesh/param_wing/VLM_mesh.msh','../mesh/param_wing/FEM_mesh.msh') 
a_new.new_wing_mesh(b,S)
vlm_mesh = '../mesh/param_wing/new_VLM.msh'
fem_mesh = '../mesh/param_wing/new_FEM.msh'
n_VLM_nodes, n_FEM_nodes, n_gamma_nodes = foff.calcule_nodes(pC)
# 3a) Computation of strains and stress
element_property,material,element_type = foff.create_dict(param_data)
my_fem = FEM_study(fem_mesh,element_type,element_property,material)
my_fem.read_mesh_file()
strain_dict, stress_dict = my_fem.get_strain_and_stress(u_est)
# 3b) u_est
my_fem.post_processing(u_est,"../results/Param_wing/u_est")
my_fem.post_processing_var1(u_est1,"../results/Param_wing/u_est1")
my_fem.post_processing_var2(u_est2,"../results/Param_wing/u_est2")
## 3c) Von Mises stress
Von_Mises_Stress = my_fem.get_Von_Mises(stress_dict)
my_fem.post_processing_Von_Mises(Von_Mises_Stress,'../results/param_wing/result_VM_iter')
print("Average Von Mises Stress = "+str(np.mean(Von_Mises_Stress))+" Pa")
stress_max = 0.
stress_max_i = 0
for i in range(len(Von_Mises_Stress)):
    stress_i = Von_Mises_Stress[i]
    if (stress_max < stress_i):
        stress_max = stress_i
        stress_max_i = i
print("Node VM max = "+str(stress_max_i))
print("Von Mises Stress max = "+str(Von_Mises_Stress[stress_max_i])+" Pa")
# 3d) g_est
my_vlm = VLM_study(vlm_mesh)
my_vlm.post_processing_gamma('Param_wing/g_est',g_est)
my_vlm.post_processing_gamma_var1('Param_wing/g_est1',g_est1)
my_vlm.post_processing_gamma_var2('Param_wing/g_est2',g_est2)
# 4) Quality of the prediction
for x in range(len(u_mean)):
    if u_mean[x] == 0:
        u_mean[x] = 1.*1e-35
tau_u = abs(np.sqrt(u_var2)/u_mean)
tau_g = abs(np.sqrt(gamma_var2)/gamma_mean)
my_fem.post_processing_erreur(tau_u,"../results/Param_wing/tau_u")
my_vlm.post_processing_erreur('Param_wing/tau_g',tau_g)
# 5) Opening GMSH files
subprocess.call(['D://gmsh-3.0.6/gmsh','../results/Param_wing/g_est.msh','../results/Param_wing/g_est1.msh','../results/Param_wing/g_est2.msh','../results/Param_wing/tau_g.msh','-open'])
subprocess.call(['D://gmsh-3.0.6/gmsh','../results/Param_wing/u_est.msh','../results/Param_wing/u_est1.msh','../results/Param_wing/u_est2.msh','../results/Param_wing/tau_u.msh','-open'])
subprocess.call(['D://gmsh-3.0.6/gmsh','../results/Param_wing/result_VM_iter.msh','-open'])