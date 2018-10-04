import numpy as np
from define_geo import airbus_wing
from create_FEM_mesh import create_FEM_mesh
from create_VLM_mesh import create_VLM_mesh
import subprocess
from PARAM_WING_VLM_FEM import PARAM_WING_VLM_FEM
from mesh_deformation_airbus import mesh_deformation_airbus
from VLM import VLM_study
from FEM import FEM_study
from transfert import transfert_matrix

"""
Definition of a function which creates the .msh files.
"""
def create_mesh_files(n_ribs_1,n_ribs_2):
    VLM_geo = '../mesh/param_wing/VLM_mesh.geo'
    a = airbus_wing()
    span_1, span_2 =  a.compute_span()
    L_1, c_1, c_2, c_3 = a.compute_cord()
    e_1,e_2,e_3 = a.compute_ep() 
    d_1 = d_2 = a.compute_diedre()
    theta_1, theta_2, theta_3 = a.compute_twist()
    phi_1, phi_2  = a.compute_phi()
    coeff = 0.70
    tck = np.array([[e_1,coeff*e_1],[e_2,coeff*e_2],[e_3,coeff*e_3]])
    chord_pos = np.array([[0.10,0.8],[0.10,0.8],[0.10,0.8]])
    create_VLM_mesh(span_1,span_2,theta_1,theta_2,theta_3,L_1,c_1,c_2,c_3,phi_1,phi_2,d_1,d_2)
    create_FEM_mesh(VLM_geo,tck,chord_pos,n_ribs_1,n_ribs_2)
    FEM_geo = '../mesh/param_wing/FEM_mesh.geo'
    #Creating a file .msh
    subprocess.call(['D://gmsh-3.0.6/gmsh', VLM_geo, '-2'])
    subprocess.call(['D://gmsh-3.0.6/gmsh', FEM_geo, '-2'])
    
"""
Definition of a function which returns the number of nodes
"""
def calcule_nodes(param_data):
    f = open('../mesh/param_wing/FEM_mesh.msh','r')
    lines = f.readlines()
    FEM_nodes = int(lines[lines.index('$Nodes\n')+1])
    f.close()
    f = open('../mesh/param_wing/VLM_mesh.msh','r')
    lines = f.readlines()
    VLM_nodes = int(lines[lines.index('$Nodes\n')+1])
    f.close()
    vlm_mesh_file = '../mesh/param_wing/new_VLM.msh'
    my_vlm = VLM_study(vlm_mesh_file,alpha = param_data[0],v_inf = param_data[1],rho = param_data[2])
    gamma_nodes = np.sum([my_vlm.surfaces[i]['points'].shape[0]for i in range(2)])
    return VLM_nodes, FEM_nodes, gamma_nodes

"""
Definition of a function which runs the solver aero_struct
"""
def run_aerostruct(param_data,b,S,phi,diedre,BF,Mach):
    #1) Applying the right values in the wing
    a_new = mesh_deformation_airbus('../mesh/param_wing/VLM_mesh.msh','../mesh/param_wing/FEM_mesh.msh') 
    a_new.new_wing_mesh(b,S,phi,diedre,BF,Mach)
    #2) Defining the files used and the parameters
    vlm_mesh_file = '../mesh/param_wing/new_VLM.msh'
    vlm_mesh_file_out = '../mesh/param_wing/new_VLM_def.msh'
    fem_mesh = '../mesh/param_wing/new_FEM.msh'
    #3) Performing the analysis and saving the results
    pw = PARAM_WING_VLM_FEM(param_data)
    U,gamma = pw.do_analysis(vlm_mesh_file,vlm_mesh_file_out,fem_mesh)
    return U, gamma

"""
Definition of a function which returns the fields needed to compute the RB-VLM parameters
"""
def VLM_params(param_data,b,S,phi,diedre,BF,Mach):
    #1) Applying the right values in the wing
    a_new = mesh_deformation_airbus('../mesh/param_wing/VLM_mesh.msh','../mesh/param_wing/FEM_mesh.msh') 
    a_new.new_wing_mesh(b,S,phi,diedre,BF,Mach)
    #2) Defining the files used and the parameters
    vlm_mesh_file = '../mesh/param_wing/new_VLM.msh'
    #3) Performing the analysis and saving the results
    vlm_part = VLM_study(vlm_mesh_file,alpha = param_data[0],v_inf = param_data[1],rho = param_data[2])
    A = vlm_part.get_A()
    B = vlm_part.get_B()
    return A, B

"""
Definition of a function which returns the Strucutural Force
"""
def get_Fs(gamma,param_data):
    vlm_mesh_file = '../mesh/param_wing/new_VLM.msh'
    fem_mesh = '../mesh/param_wing/new_FEM.msh'
    my_vlm = VLM_study(vlm_mesh_file,alpha = param_data[0],v_inf = param_data[1],rho = param_data[2])
    my_pw = PARAM_WING_VLM_FEM(param_data)
    Fa = my_vlm.compute_forces_ROB(gamma)
    vlm_nodes,fem_nodes = my_pw.get_nodes(fem_mesh,my_vlm)
    H = transfert_matrix(fem_nodes,vlm_nodes,function_type='thin_plate')
    F_a = my_vlm.compute_force_at_nodes(vlm_nodes,Fa)
    Fs = np.dot(H.T,F_a)
    return Fs

def get_FEM(param_data,F_s):
    vlm_mesh_file = '../mesh/param_wing/new_VLM.msh'
    my_vlm = VLM_study(vlm_mesh_file,alpha = param_data[0],v_inf = param_data[1],rho = param_data[2])
    element_property,material,element_type = create_dict(param_data)
    fem_mesh = '../mesh/param_wing/new_FEM.msh'
    my_fem = FEM_study(fem_mesh,element_type,element_property,material)
    elem_dict,element_tot = my_fem.read_mesh_file()
    fem_nodes_tot = elem_dict['nodes'].copy()
    elm_skins = elem_dict['element_sets'][1]['elements']
    node_skins = np.unique(elm_skins[:,3:6])
    node_skins = np.array(node_skins,dtype = int)
    fem_nodes_ind = node_skins
    my_fem.assembling_K()
    #creating the rhs for the fem study
    rhs = np.zeros((len(fem_nodes_tot)*6,))
    for i in range(3):
        rhs[6*(fem_nodes_ind-1)+i] = F_s[:,i]
    my_fem.set_rhs(rhs)    
    #Boundary condition
    my_pw = PARAM_WING_VLM_FEM(param_data)
    _,fem_nodes = my_pw.get_nodes(fem_mesh,my_vlm)
    ind = fem_nodes[:,1]<=0.06
    nodes_sets = [fem_nodes_ind[ind]]
    dof = [[0,1,2,3,4,5]]
    my_fem.boundary_conditions(nodes_sets,dof)
    K = my_fem.get_K()
    Fs = my_fem.get_rhs()
    return K, Fs

"""
Definition of a function which creates the dictionaries needed for beginning the FEM analysis
"""
def create_dict(param_data):
    element_property = {'tri':{'"skins"':{'h':param_data[5]},'"ribs"':{'h':param_data[6]},
                                   '"spars_le"':{'h':param_data[7]},'"spars_te"':{'h':param_data[8]}}}
    material = {'tri':{'"skins"':{'E':param_data[3],'nu':param_data[4]},'"ribs"':{'E':param_data[3],'nu':param_data[4]},
                       '"spars_le"':{'E':param_data[3],'nu':param_data[4]},'"spars_te"':{'E':param_data[3],'nu':param_data[4]}}}
    element_type = {'tri':'DKT'}
    return element_property,material,element_type