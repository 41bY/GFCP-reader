# -*- coding: utf-8 -*-
"""
Created on Sat Oct  1 10:22:10 2022

@author: simula
"""

import Read_input as ri
import numpy as np

###############################################################################

def test_search_input_tags():
    
    good_list = ['Line1 token tk',                  \
                 'iprint = 10',                     \
                 'dt = 4.0d0',                      \
                 'Line useless to be skipped',     \
                'celldm = 4.365354,',             \
                'nat = 5,',                       \
                'ntyp = 2',                        \
                'CELL_PARAMETERS',                \
                '4.00 0.00 0.000',                \
                '0.00 2.00 0.000',                \
                '0.00 0.00 8.000',                \
                'ATOMIC_SPECIES',                 \
                'C  12.012 c.upm',                \
                'H   2.012 h.upm',                \
                'ATOMIC_POSITIONS',               \
                'C 0.00 0.00 0.000',              \
                'C 1.00 0.00 0.000',              \
                'C 0.00 0.00 8.000',              \
                'H 1.00 0.00 0.000',              \
                'H 0.00 0.00 8.000',              \
                'Other stuff irrelevant']
        
    corrupte_tag_list = ['Line1 token tk',                  \
                 'iprint = 10',                     \
                 'dt = ',                      \
                 'Line useless to be skipped',     \
                'celldm = 4.365354,',             \
                'nat = 5,',                       \
                'ntyp = 2',                        \
                'CELL_PARAMETERS',                \
                '4.00 0.00 0.000',                \
                '0.00 2.00 0.000',                \
                '0.00 0.00 8.000',                \
                'ATOMIC_SPECIES',                 \
                'C  12.012 c.upm',                \
                'H   2.012 h.upm',                \
                'ATOMIC_POSITIONS',               \
                'C 0.00 0.00 0.000',              \
                'C 1.00 0.00 0.000',              \
                'C 0.00 0.00 8.000',              \
                'H 1.00 0.00 0.000',              \
                'H 0.00 0.00 8.000',              \
                'Other stuff irrelevant']

    no_tag_list = ['Line1 token tk',                  \
                 'iprint = 10',                     \
                 'dt = 4.0d0',                      \
                 'Line useless to be skipped',     \
                'celldm = 4.365354,',             \
                'nat = 5,',                       \
                'ntyp = 2',                        \
                'CELL_PARAMETERS',                \
                '4.00 0.00 0.000',                \
                '0.00 2.00 0.000',                \
                '0.00 0.00 8.000',                \
                'C  12.012 c.upm',                \
                'H   2.012 h.upm',                \
                'ATOMIC_POSITIONS',               \
                'C 0.00 0.00 0.000',              \
                'C 1.00 0.00 0.000',              \
                'C 0.00 0.00 8.000',              \
                'H 1.00 0.00 0.000',              \
                'H 0.00 0.00 8.000',              \
                'Other stuff irrelevant']    
    
    found = ri.search_input_tags(good_list)
    assert len(found) == 8
    assert found[0] == '10'
    assert found[1] == '4.00'
    assert found[2] == '4.365354'
    assert found[3] == '5'
    assert found[4] == '2'
    assert found[5] == [['4.00', '0.00', '0.000'], ['0.00', '2.00', '0.000'], ['0.00', '0.00', '8.000']]   
    assert found[6] == [['C','12.012'],['H', '2.012']]
    assert found[7] == ['C', 'C', 'C', 'H', 'H']
    
    found = ri.search_input_tags(corrupte_tag_list)
    assert len(found) == 8
    assert found[0] == '10'
    assert found[1] == ''
    assert found[2] == '4.365354'
    assert found[3] == '5'
    assert found[4] == '2'
    assert found[5] == [['4.00', '0.00', '0.000'], ['0.00', '2.00', '0.000'], ['0.00', '0.00', '8.000']]   
    assert found[6] == [['C','12.012'],['H', '2.012']]
    assert found[7] == ['C', 'C', 'C', 'H', 'H']
    
    found = ri.search_input_tags(no_tag_list)
    assert len(found) == 7
    assert found[0] == '10'
    assert found[1] == '4.00'
    assert found[2] == '4.365354'
    assert found[3] == '5'
    assert found[4] == '2'
    assert found[5] == [['4.00', '0.00', '0.000'], ['0.00', '2.00', '0.000'], ['0.00', '0.00', '8.000']]   
    assert found[6] == ['C', 'C', 'C', 'H', 'H'] 

###############################################################################

def test_process_tag():
    
    good_tag = 'nat = 20'
    bad_1_tag = 'nat = '
    bad_2_tag = 'nat = dcfr'
    
    assert ri.process_tag(good_tag) == '20'
    assert ri.process_tag(bad_1_tag) == ''
    assert ri.process_tag(bad_2_tag) == ''

###############################################################################

def test_process_cell():
    
    good_cell = ['CELL', '1.00 0.00 0.00', '0.00 2.00 0.00', '0.00 0.00 5.00']
    bad_1_cell = ['CELL', '1.00 0.00 0.00', '0.00 0.00', '0.00 0.00 5.00']
    bad_2_cell = ['CELL', '1.00 0.00 0.00', '0.00 2.00 0.00']
    bad_3_cell = ['CELL', '1.00 0.00 0.00', '0.00 2.00 0.00', '0.00 abcd 5.00']
    
    assert ri.process_cell(good_cell) == [['1.00', '0.00', '0.00'],['0.00', '2.00', '0.00'],['0.00', '0.00', '5.00']]
    assert ri.process_cell(bad_1_cell) == []
    assert ri.process_cell(bad_2_cell) == []
    assert ri.process_cell(bad_3_cell) == []

###############################################################################

def test_process_species():
    
    ntyp = 3
    good_specie = ['SPECIE', 'C 12.012 van.bm', 'H 2.012 van.bm', 'Fe 54.012 van.bm']
    bad_1_specie = ['SPECIE', 'C van.bm', 'H 2.012 van.bm', 'Fe 54.012 van.bm']
    bad_2_specie = ['SPECIE', 'C 12.012 van.bm', 'Fe 54.012 van.bm']
    bad_3_specie = ['SPECIE', 'C 12.012 van.bm', 'H abcde van.bm', 'Fe 54.012 van.bm']
    
    assert ri.process_species(good_specie, ntyp) == [['C', '12.012'],['H', '2.012'],['Fe', '54.012']]
    assert ri.process_species(bad_1_specie, ntyp) == []
    assert ri.process_species(bad_2_specie, ntyp) == []
    assert ri.process_species(bad_3_specie, ntyp) == []

###############################################################################

def test_process_atoms():
    
    nat = 3
    good_atoms = ['ATOMS', 'C 0 1 2', 'H 0 1 2', 'Fe 0 1 2']
    bad_1_atoms = ['ATOMS', 'C 0 1 2', 'H 0 1 2']
    bad_2_atoms = ['ATOMS', 'C 0 1 2', '', 'Fe 0 1 2']
    
    assert ri.process_atoms(good_atoms, nat) == ['C', 'H', 'Fe']
    assert ri.process_atoms(bad_1_atoms, nat) == []
    assert ri.process_atoms(bad_2_atoms, nat) == []

###############################################################################

def test_check_founds():
    
    founds = ['10', '4.0', '2.32', '2', '5', \
                    [['1.00','0.00','0.00'],['0.00','2.00','0.00'],['0.00','0.00','5.00']], \
                    [['C',' 12.012'],['H',' 2.012',],['Fe','54.012']], \
                    ['C', 'H', 'Fe']]
        
    found_list = ri.check_founds(founds)
    
    assert  isinstance(found_list[0], int)
    assert  isinstance(found_list[1], float)
    assert  isinstance(found_list[2], float)
    assert  isinstance(found_list[3], int)
    assert  isinstance(found_list[4], int)
    assert  isinstance(found_list[5], np.ndarray)
    assert  isinstance(found_list[6], list)
    assert  isinstance(found_list[7], list)