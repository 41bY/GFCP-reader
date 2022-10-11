# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 13:09:00 2022

@author: simula
"""
import numpy as np
import Read_output as ro
import io


###############################################################################

def test_check_temp():
    
    temp_good = '!     kt_gf,     nfi,          =     304.00        10'
    temp_no_tag = '!     nfi,          =     304.00        10'
    temp_corrupted = '!     kt_gf,     nfi,          ='
    
    
    assert ro.check_temp(temp_good)[0] == 304.00
    assert ro.check_temp(temp_good)[1] == ''
    assert ro.check_temp(temp_no_tag)[0] == None
    assert ro.check_temp(temp_no_tag)[1] != ''
    assert ro.check_temp(temp_corrupted)[0] == None
    assert ro.check_temp(temp_corrupted)[1] != ''
                
###############################################################################

def test_check_force():
    
    force_good = '!      lower[GPa] =    0.09767  -0.08484  -0.28874'
    force_no_tag = '!      lower[GPa]   0.09767  -0.08484  -0.28874'
    force_corrupted = '!      lower[GPa] =    -0.08484  -0.28874'
    
    
    assert ro.check_force(force_good)[0] == [0.09767, -0.08484, -0.28874]
    assert ro.check_force(force_good)[1] == ''
    assert ro.check_force(force_no_tag)[0] == None
    assert ro.check_force(force_no_tag)[1] != ''
    assert ro.check_force(force_corrupted)[0] == None
    assert ro.check_force(force_corrupted)[1] != ''     
        
###############################################################################

def test_check_qm_atoms():
    
    qm_atoms_good = ['ATOMIC_POSITIONS', 'H 1 1 1', 'C 2 2 2']
    qm_atoms_notag = ['H 1 1 1', 'C 2 2 2']
    qm_atoms_corrupted = ['ATOMIC_POSITIONS', 'H 1 1 1', 'C 2 2']

    assert ro.check_qm_atoms(qm_atoms_good)[0] == [['H', '1', '1', '1'], ['C', '2', '2', '2']]
    assert ro.check_qm_atoms(qm_atoms_good)[1] == ''
    assert ro.check_qm_atoms(qm_atoms_notag)[0] == None
    assert ro.check_qm_atoms(qm_atoms_notag)[1] != ''
    assert ro.check_qm_atoms(qm_atoms_corrupted)[0] == None
    assert ro.check_qm_atoms(qm_atoms_corrupted)[1] != '' 

###############################################################################

def test_check_gf_atoms():

    gf_atoms_good = ['GF_ATOM_POSITIONS', '1 1 1', '2 2 2']
    gf_atoms_notag = ['1 1 1', '2 2 2']
    gf_atoms_corrupted = ['GF_ATOM_POSITIONS', '1 1', '2 2 2']

    assert ro.check_gf_atoms(gf_atoms_good)[0] == [['GF', '1', '1', '1'], ['GF', '2', '2', '2']]
    assert ro.check_gf_atoms(gf_atoms_good)[1] == ''
    assert ro.check_gf_atoms(gf_atoms_notag)[0] == None
    assert ro.check_gf_atoms(gf_atoms_notag)[1] != ''
    assert ro.check_gf_atoms(gf_atoms_corrupted)[0] == None
    assert ro.check_gf_atoms(gf_atoms_corrupted)[1] != '' 

###############################################################################

# def Error_msg(iteration : int,  num_relax : int, block_size : int , num_line_block : int):
#     file_line_num = 3 + (iteration + num_relax)*block_size + num_line_block
#     file_line_num = int(file_line_num)

#     msg = ' '+str(file_line_num)+' in simulation output file.'
    
#     return msg

# ###############################################################################
                
# def initialize_output_files(path : str):
    
#     TEM_fname = path+'Temp.dat'
#     FUP_fname = path+'Forces_gf_up.dat'
#     FDW_fname = path+'Forces_gf_dw.dat'
#     POS_fname = path+'Pos.xyz'

#     #Open output files to write them in the main loop
#     f_TEM = open(TEM_fname, 'w') 
#     f_FOR_up = open(FUP_fname, 'w')
#     f_FOR_dw = open(FDW_fname, 'w')
#     f_POS = open(POS_fname, 'w')
    
#     f_TEM.writelines('Temperature of GF atoms: TimeStep(fs) T(Kelvin) \n')  
#     f_FOR_up.writelines('Total forces on upper GF atoms: TimeStep(fs) Fx Fy Fz(GPa) \n')
#     f_FOR_dw.writelines('Total forces on lower GF atoms: TimeStep(fs) Fx Fy Fz(GPa) \n')
    
#     return f_TEM, f_FOR_up, f_FOR_dw, f_POS

# ###############################################################################
                
# def get_str_values(block : list, tags_num : list, nat_qm : int, nat_gf : int):
    
#     temp_str = block[tags_num[0]]
#     fup_str = block[tags_num[1]]
#     fdw_str = block[tags_num[2]]
#     control_str = block[tags_num[3]]
#     qm_atoms_str_list = block[tags_num[4]:tags_num[4]+nat_qm+1]
#     gf_atoms_str_list = block[tags_num[5]:tags_num[5]+nat_gf+1]
    
#     return temp_str, fup_str, fdw_str, control_str, qm_atoms_str_list, gf_atoms_str_list
                   
# ###############################################################################
                
# def order_all_atoms(list_a : list, list_b : list, idx_sort : list):
    
#     list_tot = list_a + list_b
#     lis_ord = np.array(list_tot, dtype=str)
#     lis_ord = lis_ord[idx_sort]
    
#     return lis_ord

# ###############################################################################

# def get_formatted_data(time : float, temp : float, fup : list, fdw : list, atoms : np.ndarray):
#     atom_pos = atoms[:,1:].astype(float)
#     atom_label = atoms[:,0]
    
#     tempVStime = [[time, temp]]
#     fupVStime = [[time]+fup]
#     fdwVStime = [[time]+fdw]
    
#     nat = atoms[:,0].size
#     header = str(nat)+'\nATOMIC_POSITIONS (angstrom); time (fs) = %10.6f \n' %time
    
#     return tempVStime, fupVStime, fdwVStime, atom_label, atom_pos, header







