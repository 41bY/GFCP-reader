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







