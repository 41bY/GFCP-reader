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

def skip_lines(file : io.TextIOBase, lines_to_skip : int, n_line : int):
    
    #Check method functioning conditions
    if not isinstance(file, io.TextIOBase):
        raise ValueError('The object passed is not a text file object!')
        
    if not isinstance(lines_to_skip, int):
        raise ValueError('Number of lines to be skipped must be an integer!')        
    
    if not lines_to_skip > 0:
        raise ValueError('Number of lines to be skipped must be positive!')
    
    for i in range(lines_to_skip):
        line = file.readline()
        n_line += 1
    
    return line, n_line

###############################################################################

def read_tag(file : io.TextIOBase, search_list : list, n_line : int):
    
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
        n_line += 1
        
        #EOF error
        if line == '':
            line = 'Error : EOF'
            return found_dict, line, n_line
        
        for search in search_list:
            remove = None
            
            if search in line:
                tokens = line.split('=')
                
                #Empty field error
                if len(tokens) < 2:
                    line = 'Error : Empty field'
                    return found_dict, line, n_line
                
                token = tokens[1]
                found_dict[search] = token
                
                remove = search
                break
        
        #Remove found element from the search list
        if remove != None:
            search_list.remove(remove)
            
            #If search list is void return map with founds and control line
            if len(search_list) == 0: 
                return found_dict, line, n_line
            
    return found_dict, line, n_line

###############################################################################

def errors_read_tag(file_path, tags_map, line):
    
    if 'EOF' in line:
        not_found = [ele for ele in tags_map.keys() if tags_map[ele] == None]
        ErrorMSG = 'EOF reached in file ' + file_path + ' while searching for '
        for it in not_found: ErrorMSG += it+' '
        raise ValueError(ErrorMSG)
        
    elif 'Empty field' in line:
        
        for tag, ele in tags_map.items():
            
            if ele == '':
                empty = tag
                break
            
        ErrorMSG = 'In file ' + file_path + ' empty field at tag: ' + empty
        raise ValueError(ErrorMSG)
    

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
        
        line = 'Start'
        n_line = 0
        while line != '':     
            
            line = f.readline()
            
            
            if '&CONTROL' in line:
                tags = ['iprint', 'dt']
                tags_map, line = read_tag(f, tags)
                if 'Error = ' in line: errors_read_tag(file_path, tags_map, line)
                
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

def write_xyz(file : io.TextIOBase, data : tuple, header : str):
    #Check method functioning conditions
    if not isinstance(file, io.TextIOBase):
        raise ValueError('The object passed is not a text file object!')
    
    if not isinstance(data, tuple):
        raise ValueError('Data must be a tuple!')
    
    if not isinstance(header, str):
        raise ValueError('Header must be a string!')
    
    nat_pos = len(data[1][:,0])
    nat_lbl = len(data[0][:])
    if nat_pos != nat_lbl:
        raise ValueError('The row \'data\'[0] must have a size equal to the columns of \'data\'[1]!')
    
    file.writelines('%d \n' %nat_pos)
    file.writelines(header+'\n')
    for i, atom in enumerate(data[1]):
        file.writelines('%4s %25.15f %25.15f %25.15f \n' %(data[0][i], atom[0], atom[1], atom[2]))

###############################################################################                
           
def write_info(file_path : str, info_data : list):
    #Check method functioning conditions
    if not isinstance(file_path, str):
        raise ValueError('First argument must be a file path, it should be a str!')
    
    if not isinstance(info_data, list):
        raise ValueError('Info_data must be a list!')
    
    if len(info_data) != 11:
        raise ValueError('info_data must be a precise list of 11 elements!')
    
    with open(file_path, 'w') as f_info:
        
        cell = info_data[0]
        f_info.writelines('Cell (angstrom):\n')
        for i in range(3):
            f_info.writelines('%15.10f %15.10f %15.10f' %(cell[i][0], cell[i][1], cell[i][2]))
            f_info.writelines('\n')
        
        types = info_data[1]
        f_info.writelines('\nAtom types (type|mass(a.m.u.)|#):\n')
        for at in types:
            f_info.writelines('%4s %10.4f %5d' %(at[0], at[1], at[2]))
            f_info.writelines('\n')
        
        nat_tot = info_data[2]            
        f_info.writelines('\nTotal atoms number = %d \n' %nat_tot)
        
        nat_qm = info_data[3]
        f_info.writelines('\nQM atoms number = %d \n' %nat_qm)
        
        nat_gf = info_data[4]
        f_info.writelines('\nGF atoms number = %d \n' %nat_gf)
        
        data_dt = info_data[5]
        f_info.writelines('\nData time step = %10.9f fs\n' %data_dt)
        
        dt = info_data[6]
        f_info.writelines('\nSimulation time step = %10.9f fs\n' %dt)
        
        step = info_data[7]
        f_info.writelines('\nData step = %d \n' %step)
        
        data_step = info_data[8]
        f_info.writelines('\nSimulation step =  %d \n' %data_step)
        
        #NAT for down slab (the frist one) corresponds to starting index upper layer!
        f_info.writelines('\nSlabs info: (Slab \t #atoms \t Mass(a.m.u.))\n')
        
        dw_slab_num, dw_mass_slabs = info_data[9]
        f_info.writelines('Lower \t %d \t %10.5f \n' %(dw_slab_num, dw_mass_slabs))
        
        up_slab_num, up_mass_slabs = info_data[10]
        f_info.writelines('Upper \t %d \t %10.5f \n' %(up_slab_num, up_mass_slabs)) 

###############################################################################              
                    
                    
                    
                    