# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 13:09:00 2022

@author: simula
"""
import numpy as np
import Read_reference_output as rr

def test_find_tags_positions():
    
    sample_list = ['kt_gf = 300', 'irrelevant', 'upper[GPa]', 'lower[GPa]']
    
    line_num = rr.find_tags_positions(sample_list)
    
    assert line_num[0] == 0
    assert line_num[1] == 2
    assert line_num[2] == 3

###############################################################################

def test_get_qm_atoms():
    
    sample_list = ['H 1 1 1', 'C 1 3 4', 'abcd']
    
    qm_lbl, qm_pos = rr.get_qm_atoms(sample_list)
    
    assert qm_lbl == ['H', 'C']
    assert qm_pos == [['1', '1', '1'], ['1', '3', '4']]

###############################################################################

def test_get_gf_atoms():
    
    sample_list = ['1 1 1', '1 3 4', 'abcd']
    
    num, gf_pos = rr.get_gf_atoms(sample_list)
    
    assert num == len(gf_pos)
    assert num == 2
    assert gf_pos == [['1','1','1'],['1','3','4']]
          
###############################################################################

def test_sort_atoms():
    
    list2sort = [['1','2','4'], ['0','1','1'], ['3','2','-5']]
    
    idx_z_sort, up_slab_index = rr.sort_atoms(list2sort)
    
    assert (idx_z_sort == np.array([2, 1, 0], dtype=int)).all
    assert up_slab_index == 1
