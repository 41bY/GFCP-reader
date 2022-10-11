# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 13:09:00 2022

@author: simula
"""
import numpy as np
import Read_reference_output as rr

# def read_ref_output(file_path : str):
    
#     file = open(file_path, 'r')
    
#     #Skip comment lines
#     for i in range(3): next(file)
    
#     #Check starting iteration
#     ref_tag = 'kt_gf'
#     line = next(file)
#     if not (ref_tag in line): raise ValueError('Output file format is not as expected!')
    
#     #Read an iteration block
#     block = read_ref_block(file, ref_tag)
#     block = [line] + block
#     block_size = len(block)
#     file.close()
            
#     #Get the line numbers of the relevant tags: 0->Temp, 1->Fup, 2->Fdw, 3->after, 4->qm, 5->gf
#     tags_line_num = find_tags_positions(block)
#     if len(tags_line_num) != 6: raise ValueError('Output file format is not as expected!')
    
#     #Get the atomic label
#     qm_block = block[tags_line_num[4]+1:]
#     qm_label, qm_pos = get_qm_atoms(qm_block)
    
#     #Get gf atom number
#     gf_block = block[tags_line_num[5]+1:]
#     nat_gf, gf_pos = get_gf_atoms(gf_block)
    
#     #Sort all atoms along z
#     at_pos = qm_pos + gf_pos
#     idx_z_sort, up_slab_idx = sort_atoms(at_pos)
    
#     return tags_line_num, block_size, qm_label, nat_gf, idx_z_sort, up_slab_idx

###############################################################################

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
