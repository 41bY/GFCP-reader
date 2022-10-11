# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 13:09:00 2022

@author: simula
"""
import numpy as np
import io

def read_ref_output(file_path : str):
    """
    'main' program of this module. It reads a 'gfcp.out' output file and it controls
    the other methods below. Open the 'gfcp.out' output file and read it in block
    which is controlled by the tag 'kt_gf'. This is done because 'gfcp.out' is
    usually a very large file but it is structured in 'iteration block' and each
    block has the same format of the others.
    The objective of this method is to read a reference block in order to get
    the indeces of each relevant tag within the block and thus speeding up the
    'real' reading.

    Parameters
    ----------
    file_path : str
        'gfcp.out' output file path.

    Raises
    ------
    errors
        Inside 'check_founds' if the 'cp.x' is corrupted.

    Returns
    -------
    found_list : list
        Output list containing the values associated to each tag of the proper
        type depending on the tag: 
        'iprint' -> int; 'dt' -> float; 'celldm' -> float; 'nat' -> int; 'ntyp' -> int;
        'CELL_PARAMETERS' -> np.array(float); 'ATOMIC_SPECIES' -> list of [str, float];
        'ATOMIC_POSITIONS' -> list of str

    Parameters
    ----------
    file_path : str
        DESCRIPTION.

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    tags_line_num : TYPE
        DESCRIPTION.
    block_size : TYPE
        DESCRIPTION.
    qm_label : TYPE
        DESCRIPTION.
    nat_gf : TYPE
        DESCRIPTION.
    idx_z_sort : TYPE
        DESCRIPTION.
    up_slab_idx : TYPE
        DESCRIPTION.

    """
    
    #Open file
    file = open(file_path, 'r')
    
    #Skip comment lines
    for i in range(3): next(file)
    
    #Check starting iteration
    ref_tag = 'kt_gf'
    line = next(file)
    if not (ref_tag in line): raise ValueError('Output file format is not as expected!')
    
    #Read an iteration block
    block = read_ref_block(file, ref_tag)
    block = [line] + block
    block_size = len(block)
    file.close()
            
    #Get the line numbers of the relevant tags: 0->Temp, 1->Fup, 2->Fdw, 3->after, 4->qm, 5->gf
    tags_line_num = find_tags_positions(block)
    if len(tags_line_num) != 6: raise ValueError('Output file format is not as expected!')
    
    #Get the atomic label
    qm_block = block[tags_line_num[4]+1:]
    qm_label, qm_pos = get_qm_atoms(qm_block)
    
    #Get gf atom number
    gf_block = block[tags_line_num[5]+1:]
    nat_gf, gf_pos = get_gf_atoms(gf_block)
    
    #Sort all atoms along z
    at_pos = qm_pos + gf_pos
    idx_z_sort, up_slab_idx = sort_atoms(at_pos)
    
    return tags_line_num, block_size, qm_label, nat_gf, idx_z_sort, up_slab_idx

###############################################################################

def read_ref_block(file : io.TextIOBase, ref_tag : str):
    
    block = []
    for line in file:
        
        if ref_tag in line: break
        else: block.append(line)
    
    return block

###############################################################################

def find_tags_positions(input_list : list):
    
    line_num = []
    for n, line in enumerate(input_list):
        
        if 'kt_gf' in line: 
            line_num.append(n)
        
        elif 'upper[GPa]' in line: 
            line_num.append(n)            
        
        elif 'lower[GPa]' in line: 
            line_num.append(n)      
        
        elif 'after' in line: 
            line_num.append(n)      
            
        elif 'ATOMIC_POSITIONS' in line: 
            line_num.append(n)

        elif 'GF_ATOM_POSITIONS' in line: 
            line_num.append(n)
    
    return line_num

###############################################################################

def get_qm_atoms(input_list : list):
    
    qm_lbl = []
    qm_pos = []
    for line in input_list:
        tokens = line.split()
        
        if len(tokens) != 4: break
        
        else:
            label = tokens[0]
            cord = tokens[1:]
            qm_lbl.append(label)
            qm_pos.append(cord)
    
    return qm_lbl, qm_pos

###############################################################################

def get_gf_atoms(input_list : list):
    
    gf_pos = []
    for n, line in enumerate(input_list):
        tokens = line.split()
        
        if len(tokens) != 3: break
        gf_pos.append(tokens)
    
    return n, gf_pos
          
###############################################################################

def sort_atoms(input_list : list):
    
    z_pos = [atom[2] for atom in input_list]
    z_pos = np.array(z_pos, dtype=float)
    
    #Sort along z
    idx_z_sort = np.argsort(z_pos)
    sorted_z_pos = z_pos[idx_z_sort]
    
    #Get index maximum z separation
    dz = sorted_z_pos[1:] - sorted_z_pos[:-1]
    up_slab_index = np.argmax(dz) + 1
    
    return idx_z_sort, up_slab_index      