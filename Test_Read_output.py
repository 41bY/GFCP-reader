# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 13:09:00 2022

@author: simula
"""
import Read_output as ro

###############################################################################

def test_check_temp_good():
    """
    Test the success of the function to extract the temperature data (float)
    from the temperature line str.
    0 index is for the value
    1 index is for the error message
    

    Returns
    -------
    None.

    """
    
    temp_line = '!     kt_gf,     nfi,          =     304.00        10'
    
    assert ro.check_temp(temp_line)[0] == 304.00
    assert ro.check_temp(temp_line)[1] == ''

###############################################################################

def test_check_temp_notag():
    """
    Test the success of the function to recognize the missing tag from the
    temperature line str.
    0 index is for the value
    1 index is for the error message    

    Returns
    -------
    None.

    """
    
    temp_line = '!     nfi,          =     304.00        10'
    
    assert ro.check_temp(temp_line)[0] == None
    assert ro.check_temp(temp_line)[1] == '\'kt_gf\' tag corrupted at line'
    
###############################################################################

def test_check_temp_corrupted():
    """
    Test the success of the function to recognize the corrupted field in the
    temperature line str.
    0 index is for the value
    1 index is for the error message    

    Returns
    -------
    None.

    """
    
    temp_line = '!     kt_gf,     nfi,          ='
    
    assert ro.check_temp(temp_line)[0] == None
    assert ro.check_temp(temp_line)[1] == '\'kt_gf\' tag corrupted at line'
                
###############################################################################

def test_check_force_good():
    """
    Test the success of the function to extract the force data (list of float)
    from the force line str
    0 index is for the value
    1 index is for the error message    

    Returns
    -------
    None.

    """
    
    force_line = '!      lower[GPa] =    0.09767  -0.08484  -0.28874'
    
    
    assert ro.check_force(force_line)[0] == [0.09767, -0.08484, -0.28874]
    assert ro.check_force(force_line)[1] == ''

###############################################################################

def test_check_force_notag():
    """
    Test the success of the function to recognize the corrupted field ('=' missing)
    in the force line str.
    0 index is for the value
    1 index is for the error message    

    Returns
    -------
    None.

    """    
    
    force_line = '!      lower[GPa]   0.09767  -0.08484  -0.28874'
    
    assert ro.check_force(force_line)[0] == None
    assert ro.check_force(force_line)[1] == '\'Force\' tag corrupted at line'

###############################################################################

def test_check_force_corrupted():
    """
    Test the success of the function to recognize the corrupted field (only 2 fields)
    in the force line str.
    0 index is for the value
    1 index is for the error message    

    Returns
    -------
    None.

    """
    
    force_line = '!      lower[GPa] =    -0.08484  -0.28874'
    
    assert ro.check_force(force_line)[0] == None
    assert ro.check_force(force_line)[1] == '\'Force\' tag corrupted at line'
        
###############################################################################

def test_check_qm_atoms_good():
    """
    Test the success of the function to extract the qm atom data (list of list of str)
    from the qm atom list.
    0 index is for the value
    1 index is for the error message    

    Returns
    -------
    None.

    """
    
    qm_atoms_list = ['ATOMIC_POSITIONS', 'H 1 1 1', 'C 2 2 2']

    assert ro.check_qm_atoms(qm_atoms_list)[0] == [['H', '1', '1', '1'], ['C', '2', '2', '2']]
    assert ro.check_qm_atoms(qm_atoms_list)[1] == ''

###############################################################################

def test_check_qm_atoms_notag():
    """
    Test the success of the function to recognize the missing tag from the
    qm atom list.
    0 index is for the value
    1 index is for the error message    

    Returns
    -------
    None.

    """
    
    qm_atoms_list = ['H 1 1 1', 'C 2 2 2']

    assert ro.check_qm_atoms(qm_atoms_list)[0] == None
    assert ro.check_qm_atoms(qm_atoms_list)[1] == '\'ATOMIC_POSITIONS\' tag corrupted at line'
    
###############################################################################

def test_check_qm_atoms_corrupted():
    """
    Test the success of the function to recognize the corrupted field from the
    qm atom list.
    0 index is for the value
    1 index is for the error message    

    Returns
    -------
    None.

    """
    
    qm_atoms_list = ['ATOMIC_POSITIONS', 'H 1 1 1', 'C 2 2']

    assert ro.check_qm_atoms(qm_atoms_list)[0] == None
    assert ro.check_qm_atoms(qm_atoms_list)[1] == '\'ATOMIC_POSITIONS\' tag corrupted at line' 

###############################################################################

def test_check_gf_atoms_good():
    """
    Test the success of the function to extract the gf atom data (list of list of str)
    from the gf atom list
    0 index is for the value
    1 index is for the error message    

    Returns
    -------
    None.

    """   

    gf_atoms_list = ['GF_ATOM_POSITIONS', '1 1 1', '2 2 2']

    assert ro.check_gf_atoms(gf_atoms_list)[0] == [['GF', '1', '1', '1'], ['GF', '2', '2', '2']]
    assert ro.check_gf_atoms(gf_atoms_list)[1] == ''

###############################################################################

def test_check_gf_atoms_notag():
    """
    Test the success of the function to recognize the missing tag from the
    gf atom list.
    0 index is for the value
    1 index is for the error message    

    Returns
    -------
    None.

    """  

    gf_atoms_list = ['1 1 1', '2 2 2']

    assert ro.check_gf_atoms(gf_atoms_list)[0] == None
    assert ro.check_gf_atoms(gf_atoms_list)[1] == '\'GF_ATOM_POSITIONS\' tag corrupted at line'
    
###############################################################################

def test_check_gf_atoms_corrupted():
    """
    Test the success of the function to recognize the corrupted field from the
    gf atom list.
    0 index is for the value
    1 index is for the error message    

    Returns
    -------
    None.

    """    

    gf_atoms_list = ['GF_ATOM_POSITIONS', '1 1', '2 2 2']

    assert ro.check_gf_atoms(gf_atoms_list)[0] == None
    assert ro.check_gf_atoms(gf_atoms_list)[1] == '\'GF_ATOM_POSITIONS\' tag corrupted at line'
    
###############################################################################

def test_Error_msg_line():
    """
    Test the success of the function to assign the correct line number where the
    error occurred.

    Returns
    -------
    None.

    """
    
    err_msg = ro.Error_msg_line(iteration = 5,  num_relax = 1, block_size = 100 , \
                           num_line_block = 3)
    
    assert '606' in err_msg

###############################################################################

def test_get_str_values():
    """
    Test the success of the function to extract the lines where the relevant
    tags are (supposed to be) found based on the indeces given.

    Returns
    -------
    None.

    """
    
    block = ['Temperature', \
             'irrelevant',  \
             'irrelevant',  \
             'Force_up',    \
             'Force_dw',    \
             'irrelevant',  \
             'irrelevant',  \
             'control',     \
             'atom_pos',    \
             'H 1 1 1 ',    \
             'C 2 2 2 ',    \
             'gf_pos  ',    \
             ' 0 0 0  ',    \
             ]
    
    tags_num = [0, 3, 4, 7, 8, 11]
    
    lines = ro.get_str_values(block, tags_num, nat_qm = 2, nat_gf = 1)
    
    assert lines[0] == 'Temperature'
    assert lines[1] == 'Force_up'
    assert lines[2] == 'Force_dw'
    assert lines[3] == 'control'
    assert lines[4] == ['atom_pos', 'H 1 1 1 ', 'C 2 2 2 ']
    assert lines[5] == ['gf_pos  ', ' 0 0 0  ']

###############################################################################

def test_order_all_atoms():
    """
    Test the success of the function to merge two lists and order them based on
    a given list of indeces.

    Returns
    -------
    None.

    """
    
    list_a = ['2', '3', '4']
    list_b = ['9', '1', '8']
    idx_sort = [4, 0, 1, 2, 5, 3]
    
    ordered_list = ro.order_all_atoms(list_a, list_b, idx_sort)
    
    assert (ordered_list == ['1','2','3','4','8','9']).all

###############################################################################
