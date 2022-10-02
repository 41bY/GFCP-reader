# -*- coding: utf-8 -*-
"""
Created on Sat Oct  1 10:22:10 2022

@author: simula
"""
import numpy as np
import os
import io
from units import PhysConstants as phy

###############################################################################

def trim_string(in_str : str, mode = 'num'):
    """
    Trim a string depending on 'mode' values. Default 'mode' is set to 'num'
    which stands for numeric. Other default implementation is 'alpha' which
    stands for alphabetic. Eventually one can give a user define mode by specifing
    what to include in the trimming process.

    Parameters
    ----------
    in_str : str
        Input string to be trimmed.
    mode : str, optional
        Trimming mode. The default is 'num' and stands for numeric trimming
        where numbers and '.' are conserved in the string. The other default
        implementation is 'alpha' which stands for alphabetic where only letters
        are conserved in the string. By passing a generic string to mode the
        trimming will be personalized based on that string.

    Returns
    -------
    out_str : str
        The trimmed string starting by 'in_str' based on 'mode'.

    """
    
    if mode == 'num':
        out_str = ''.join([c for c in in_str if c in '0123456789.'])
    
    elif mode == 'alpha':
        out_str = ''.join([c for c in in_str if c in 'abcdefghkjilmnopqrstuvwxyz'])
    
    else:
        out_str = ''.join([c for c in in_str if c in mode])
    
    return out_str

###############################################################################

def check_file_existance(file_path : str):
    """
    Check the existance of a file.

    Parameters
    ----------
    file_path : str
        Name of the file and corresponding path using bash convention i.e. 
        use '/' to separate various directories and '\' for special characters.

    Raises
    ------
    ValueError
        If the file is not found.

    Returns
    -------
    bool
        True if the file is found.

    """
    if not os.path.isfile(file_path):
        raise ValueError('%s not found' %file_path)
        
    else:
        return True
    
###############################################################################

def create_dir(dir_path : str):
    """
    Check the existance of a directory, create the directory and all the 
    directories above it, if not present.

    Parameters
    ----------
    dir_path : str
        Path to the directory using bash convention i.e. use '/' to separate 
        various directories and '\' for special characters.

    Returns
    -------
    dir_path : str
        Path to the wanted directory/ies.

    """
    if not os.path.isdir(dir_path):
    
        path_to_probe = ''
        for token in dir_path.split('/'):
            path_to_probe = path_to_probe + token + '/'
            
            if not os.path.isdir(path_to_probe):
                os.mkdir(path_to_probe)
    
    return dir_path

###############################################################################

def skip_lines(file : io.TextIOBase, lines_to_skip : int):
    
    #Check method functioning conditions
    if not isinstance(file, io.TextIOBase):
        raise ValueError('The object passed is not a text file object!')
        
    if not isinstance(lines_to_skip, int):
        raise ValueError('Number of lines to be skipped must be an integer!')        
    
    if not lines_to_skip > 0:
        raise ValueError('Number of lines to be skipped must be positive!')
        
    for i in range(lines_to_skip):
        line = file.readline()
    
    return line

###############################################################################

def read_tag(file : io.TextIOBase, search_list : list):
    
    #Check method functioning conditions
    if not isinstance(file, io.TextIOBase):
        raise ValueError('The object passed is not a text file object!')
    
    if len(search_list) == 0:
        raise ValueError('The search list cannot be empty at the beginning!')
    
    for search in search_list:
        if not isinstance(search, str):
            raise ValueError('The search list must contain only strings!')
        
    
    found_dict = dict.fromkeys(search_list)
    
    line = 'Start'
    while line != '':

        line = file.readline()
        if line == '': return found_dict, line
        
        for search in search_list:
            remove = None
            
            if search in line:
                token = line.split('=')[1]
                found_dict[search] = token
                
                remove = search
                break
        
        if remove != None: #Remove found elements from the search list
            search_list.remove(remove)        
            if len(search_list) == 0: #If search list is void return map with founds
                return found_dict, line
            
    return found_dict, line

###############################################################################

def read_list(file : io.TextIOBase, num_line : int = 0, check_consistency = True):
    
    #Check method functioning conditions
    if not isinstance(file, io.TextIOBase):
        raise ValueError('The object passed is not a text file object!')
    
    if not isinstance(num_line, int):
        raise ValueError('\'num_line\' parameter must be integer!')
    
    if (num_line == 0) & (not check_consistency):
        raise ValueError('\'num_line\' parameter must be greater than zero if \
                         consistency reading is disabled!')
    
    
    out_list = []
    
    if num_line > 1: #Read 'num_line' lines
        
        begin = True
        for n in range(num_line):
            
            line = file.readline()
            if line == '': return out_list, line
            
            token = line.split()
            
            if check_consistency:
                if begin: 
                    ref = len(token)
                    begin = False
                else:
                    n_token = len(token)
                    if n_token != ref:
                        raise ValueError('Inconsistent number of columns in read list: \
                                         check the read file or the number of lines to\
                                         be read')                     
            
            out_list.append(token)
        
        return out_list, line
            
            
    elif num_line == 0: #Read lines consistently
        
        begin = True
        line = 'Start'
        while line != '':
            
            line = file.readline()
            if line == '': return out_list, line
            
            token = line.split()

            if begin: 
                ref = len(token)
                begin = False
            else:
                n_token = len(token)
                if n_token != ref:
                    break
                                          
            out_list.append(token)
        
        return out_list, line

###############################################################################

def read_CP_input(file_path : str):
    
    with open(file_path, 'r') as f:
        for line in f:     
            if line == '': break
            
            
            if '&CONTROL' in line:
                tags = ['iprint', 'dt']
                tags_map, line = read_tag(f, tags)
                if line == '': 
                    raise ValueError('Simulation input file is missing: ' %tags)
                
                iprint = tags_map['iprint']
                iprint = trim_string(iprint, mode='num')
                iprint = int(iprint)
                
                dt = tags_map['dt']
                dt = trim_string(dt, mode='num')
                dt = float(dt)*2.4189*0.01 #fs
                    
            
            elif '&SYSTEM' in line:
                tags = ['celldm', 'nat', 'ntyp']
                tags_map, line = read_tag(f, tags)
                if line == '': 
                    raise ValueError('Simulation input file is missing: ' %tags)
                
                celldm = tags_map['celldm']
                celldm = trim_string(celldm, mode='num')
                celldm = float(celldm)
                

                nat = tags_map['nat']
                nat = trim_string(nat, mode='num')
                nat = int(nat)
                
                ntyp = tags_map['ntyp']
                ntyp = trim_string(ntyp, mode='num')
                ntyp = int(ntyp)
        
        
            elif ('CELL_PARAMETERS' in line):
                token = line.split()
                
                #Read the unit length and convert it in Angstrom
                if len(token) > 1: #User defined unit length
                    unit = token[1]
                    unit = trim_string(unit, mode='alpha')
                    
                    if unit == 'bohr':
                        cell_unit = phy.bohr_to_A
                    
                    elif unit == 'alat':
                        cell_unit = celldm*phy.bohr_to_A
                    
                    elif unit == 'angstrom':
                        cell_unit = 1.0
                        
                    else:
                        raise ValueError('Improper unit used in CELL_PARAMETERS\
                                         : %s' %unit)
                        
                else: #QE default unit length
                    cell_unit = celldm*phy.bohr_to_A
                
                #Read the cell parameters              
                cell, line = read_list(f, num_line = 3, check_consistency = True)
                cell = cell_unit*np.array(cell, dtype=float)
                
                
            elif ('ATOMIC_SPECIES' in line):
                atom_species, line = read_list(f, num_line = ntyp, check_consistency = True)
                
                #Create atom types list: [kind, mass, number]
                atom_kinds = [(at[0], [at[0], float(at[1]), 0]) for at in atom_species]
                types = dict(atom_kinds)
    
                
            elif ('ATOMIC_POSITIONS' in line):
                atom_pos, line = read_list(f, num_line = nat, check_consistency = False)
                
                #Counts the number of each different atom and insert its number
                for at in atom_pos:
                    types[at[0]][2] += 1                  
                
                # types = list(types.values())

    return iprint, dt, celldm, nat, ntyp, cell, types                           
                    
###############################################################################



###############################################################################                
                    
                    
                    
                    
                    