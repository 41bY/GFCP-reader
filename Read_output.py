# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 13:09:00 2022

@author: simula
"""
import numpy as np
import io

###############################################################################

def read_block(file : io.TextIOBase, num_lines : int):
    """
    Read iteration block from GFCP output file consisting in 'num_lines'.

    Parameters
    ----------
    file : io.TextIOBase
        File object for GFCP output file
    num_lines : int
        Number of lines to be taken for a block.

    Returns
    -------
    block : list of str
        List containing the lines of the read block. If EOF is reached the lines
        are read with ''.

    """
    
    block = []
    for i in range(num_lines):
        line = file.readline()
        block.append(line)
    
    return block

###############################################################################

def check_temp(temp_str : str):
    """
    Check the line associated to temperature which contains the tag 'kt_gf'.
    Control if the line is corrupted, in that case return an error_msg.

    Parameters
    ----------
    temp_str : str
        Line associated to temperature

    Returns
    -------
    temp : float or None
        Temperature or None if line is corrupted
    err : str
        '' if no corruption occurs, error message otherwise.

    """
    
    err_msg = '\'kt_gf\' tag corrupted at line'
    
    if not ('kt_gf' in temp_str):
        err = err_msg
        temp = None
        return temp, err
    
    
    sides = temp_str.split('=')
    if len(sides) != 2: 
        err = err_msg
        temp = None
        return temp, err
    
    else:
        right_hand = sides[1].split()
        
        if len(right_hand) != 2:
            err = err_msg
            temp = None
            return temp, err
        else:
            temp = float(right_hand[0])
            err = ''
            return temp, err
                
###############################################################################

def check_force(force_str : str):
    """
    Check the line associated to forces which contains the tag 'upper[GPa]' or
    'lower[GPa]'.
    Control if the line is corrupted, in that case return an error_msg.

    Parameters
    ----------
    force_str : str
        Line associated to force

    Returns
    -------
    force : list of float or None
        Forces or None if line is corrupted
    err : str
        '' if no corruption occurs, error message otherwise.

    """
    
    err_msg = '\'Force\' tag corrupted at line'
    
    if (not ('upper[GPa]' in force_str)) & (not ('lower[GPa]' in force_str)):
        err = err_msg
        force = None
        return force, err
    
    sides = force_str.split('=')
    if len(sides) != 2: 
        err = err_msg
        force = None
        return force, err
    
    else:
        right_hand = sides[1].split()
        
        if len(right_hand) != 3:
            err = err_msg
            force = None
            return force, err
        else:
            force = [float(tk) for tk in right_hand]
            err = ''
            return force, err        
        
###############################################################################

def check_qm_atoms(atom_list : list):
    """
    Check the list associated to qm_atoms which contains the tag 'ATOMIC_POSITIONS'.
    Control if the lines are corrupted, in that case return an error_msg.

    Parameters
    ----------
    atom_list : list of str
        List associated to atoms

    Returns
    -------
    qm_atoms : list of list of str or Null list
        Lines of atom labels and coordinates or Null if corruption occurs
    err : str
        '' if no corruption occurs, error message otherwise.

    """

    err_msg = '\'ATOMIC_POSITIONS\' tag corrupted at line'
    
    if not ('ATOMIC_POSITIONS' in atom_list[0]):
        err = err_msg
        line = 0
        qm_atoms = None
        return qm_atoms, err, line
    
    else:
        qm_atoms = []
        for line, row in enumerate(atom_list[1:]):
            atom = row.split()
            
            if len(atom) != 4:
                err = err_msg
                line = line + 1
                qm_atoms = None
                return qm_atoms, err, line
            
            else:
                qm_atoms.append(atom)
        
        err = ''
        line = line + 1
        return qm_atoms, err, line

###############################################################################

def check_gf_atoms(atom_list : list):
    """
    Check the list associated to gf_atoms which contains the tag 'GF_ATOM_POSITIONS'.
    Control if the lines are corrupted, in that case return an error_msg.

    Parameters
    ----------
    atom_list : list of str
        List associated to atoms

    Returns
    -------
    gf_atoms : list of list of str or Null list
        Lines of atom labels and coordinates or Null if corruption occurs
    err : str
        '' if no corruption occurs, error message otherwise.

    """

    err_msg = '\'GF_ATOM_POSITIONS\' tag corrupted at line'
    
    if not ('GF_ATOM_POSITIONS' in atom_list[0]):
        err = err_msg
        line = 0
        gf_atoms = None
        return gf_atoms, err, line
    
    else:
        gf_atoms = []
        for line, row in enumerate(atom_list[1:]):
            atom = row.split()
            
            if len(atom) != 3:
                err = err_msg
                line = line + 1
                gf_atoms = None
                return gf_atoms, err, line
            
            else:
                atom = ['GF'] + atom
                gf_atoms.append(atom)
        
        err = ''
        line = line + 1
        return gf_atoms, err, line

###############################################################################

def write_xyz(file : io.TextIOBase, data_label : list, data_xyz : list, header : str):
    """
    Write a standard .xyz file format.

    Parameters
    ----------
    file : io.TextIOBase
        File object to be written.
    data_label : list of str
        Labels of the atomic data.
    data_xyz : list of float
        Coordinates of the atomic data.
    header : str
        Header before each iteration.

    Returns
    -------
    None.

    """
    
    file.write(header)
    for i, xyz in enumerate(data_xyz):
        file.write('%3s %25.12f %25.12f %25.12f \n' %(data_label[i], xyz[0], xyz[1], xyz[2]))

###############################################################################

def Error_msg_line(iteration : int,  num_relax : int, block_size : int , num_line_block : int):
    """
    Construct an error message indicating the line where the error has occured.

    Parameters
    ----------
    iteration : int
        Number of 'writing' iterations
    num_relax : int
        Number of skipped iterations
    block_size : int
        Number of line in a block
    num_line_block : int
        Line index inside the block where the error occurred.

    Returns
    -------
    msg : str
        Error message specifing where the error has occurred.

    """
    file_line_num = 3 + (iteration + num_relax)*block_size + num_line_block
    file_line_num = int(file_line_num)

    msg = ' '+str(file_line_num)+' in simulation output file.'
    
    return msg

###############################################################################
                
def initialize_output_files(path : str):
    """
    Given a 'path' open the different files that has to be written.
    Write header lines on them.

    Parameters
    ----------
    path : str
        Path of the output files.

    Returns
    -------
    f_TEM : io.TextIOBase
        Temperature file object.
    f_FOR_up : io.TextIOBase
        Forces up file object.
    f_FOR_dw : io.TextIOBase
        Forces down file object.
    f_POS : io.TextIOBase
        xyz file object.

    """
    
    TEM_fname = path+'Temp.dat'
    FUP_fname = path+'Forces_gf_up.dat'
    FDW_fname = path+'Forces_gf_dw.dat'
    POS_fname = path+'Pos.xyz'

    #Open output files to write them in the main loop
    f_TEM = open(TEM_fname, 'w') 
    f_FOR_up = open(FUP_fname, 'w')
    f_FOR_dw = open(FDW_fname, 'w')
    f_POS = open(POS_fname, 'w')
    
    f_TEM.writelines('Temperature of GF atoms: TimeStep(fs) T(Kelvin) \n')  
    f_FOR_up.writelines('Total forces on upper GF atoms: TimeStep(fs) Fx Fy Fz(GPa) \n')
    f_FOR_dw.writelines('Total forces on lower GF atoms: TimeStep(fs) Fx Fy Fz(GPa) \n')
    
    return f_TEM, f_FOR_up, f_FOR_dw, f_POS

###############################################################################
                
def get_str_values(block : list, tags_num : list, nat_qm : int, nat_gf : int):
    """
    Given the indeces associated to each tags, slice the block_list accordingly
    and assign the slices to each variables.

    Parameters
    ----------
    block : list of str
        Block list containing the lines of the block iteration.
    tags_num : list of int
        Indeces within the block associated to each tag.
    nat_qm : int
        Number of lines to slice when the tag is 'ATOMIC_POSITIONS'
    nat_gf : int
        Number of lines to slice when the tag is 'GF_ATOM_POSITIONS'

    Returns
    -------
    temp_str : str
        Line containing temperature.
    fup_str : str
        Line containing Fup.
    fdw_str : str
        Line containing Fdw.
    control_str : str
        Line containing control string.
    qm_atoms_str_list : list of str
        List containing qm_atoms.
    gf_atoms_str_list : list of str
        List containing gf_atoms.

    """
    
    temp_str = block[tags_num[0]]
    fup_str = block[tags_num[1]]
    fdw_str = block[tags_num[2]]
    control_str = block[tags_num[3]]
    qm_atoms_str_list = block[tags_num[4]:tags_num[4]+nat_qm+1]
    gf_atoms_str_list = block[tags_num[5]:tags_num[5]+nat_gf+1]
    
    return temp_str, fup_str, fdw_str, control_str, qm_atoms_str_list, gf_atoms_str_list
                   
###############################################################################
                
def order_all_atoms(list_a : list, list_b : list, idx_sort : list):
    """
    Join qm and gf atoms in a unique list and order them.

    Parameters
    ----------
    list_a : list of str
        List of qm atoms
    list_b : list of str
        List of gf atoms
    idx_sort : np.ndarray of int
        List of indeces to sort atoms as the reference one.

    Returns
    -------
    lis_ord : list of str
        List of all the atoms (qm+gf) ordered.

    """
    
    list_tot = list_a + list_b
    lis_ord = np.array(list_tot, dtype=str)
    lis_ord = lis_ord[idx_sort]
    
    return lis_ord

###############################################################################

def get_formatted_data(time : float, temp : float, fup : list, fdw : list, atoms : np.ndarray):
    """
    Format the data types to be written to file. Join lists of Temperature and
    Forces data with time data and prepare header for writing .xyz.

    Parameters
    ----------
    time : float
        Data step time
    temp : float
        Temperature at that step
    fup : list
        Force up at that step
    fdw : list
        DESCRIPTION.
    atoms : np.ndarray
        DESCRIPTION.

    Returns
    -------
    tempVStime : List of list of float
        1x2 matrix containing time and temperature
    fupVStime : List of list of float
        1x2 matrix containing time and Force up
    fdwVStime : List of list of float
        1x2 matrix containing time and Force down
    atom_label : np.ndarray of str
        Containing labels of each atom
    atom_pos : np.ndarray of float
        natx3 matrix containing x,y,z coordinates of each atom
    header : str
        Containing header for the .xyz file

    """
    
    tempVStime = [[time, temp]]
    fupVStime = [[time]+fup]
    fdwVStime = [[time]+fdw]
    
    atom_pos = atoms[:,1:].astype(float)
    atom_label = atoms[:,0]
    
    nat = atom_label.size
    header = str(nat)+'\nATOMIC_POSITIONS (angstrom); time (fs) = %10.6f \n' %time
    
    return tempVStime, fupVStime, fdwVStime, atom_label, atom_pos, header
