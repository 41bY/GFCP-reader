# -*- coding: utf-8 -*-
"""
Created on Sat Oct  1 10:22:10 2022

@author: simula
"""
import os
import io

###############################################################################
#Conversion factor

hatime2fs = 0.024189 #Hartree time to femtoseconds
bohr2angstrom = 0.529177 #Bohr to angstrom

###############################################################################

def trim_string(in_str : str, mode = 'alpha-num', trim_from = ''):
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
    if mode == 'alpha-num':
        out_str = ''.join([c for c in in_str if c in '0123456789.abcdefghkjilmnopqrstuvwxyz'])    
    
    elif mode == 'num':
        out_str = ''.join([c for c in in_str if c in '0123456789.+-'])
    
    elif mode == 'alpha':
        out_str = ''.join([c for c in in_str if c in 'abcdefghkjilmnopqrstuvwxyz'])
    
    else:
        out_str = ''.join([c for c in in_str if c in trim_from])
    
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
    n_line = 0
    while line != '':

        line = file.readline()
        n_line += 1
        
        #EOF error
        if line == '':
            line = 'ERROR_read_tag = EOF'
            return found_dict, line, n_line
        
        for search in search_list:
            remove = None
            
            if search in line:
                tokens = line.split('=')
                
                #Corrupted field error
                if len(tokens) < 2:
                    line = 'ERROR_read_tag = Corrupted field'
                    found_dict[search] = ''
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

def errors_read_tag(file_path, tags_map, line, n_line):
    
    if 'EOF' in line:
        not_found = [ele for ele in tags_map.keys() if tags_map[ele] == None]
        not_found_MSG = ''
        for nf in not_found: not_found_MSG += nf+' '
        
        ErrorMSG = 'EOF reached in file \'' + file_path + '\' at line '+ str(n_line) \
            + ': ' + not_found_MSG + 'not found!'
        raise ValueError(ErrorMSG)
        
    elif 'Corrupted field' in line:
        for tag, ele in tags_map.items():
            if ele == '':
                empty = tag
                break
            
        ErrorMSG = 'In file \'' + file_path + '\' ,corrupted field for tag: \'' + empty + \
            '\' at line ' + str(n_line)
        raise ValueError(ErrorMSG)
    

###############################################################################

def read_list(file : io.TextIOBase, num_line : int = 0, num_col : int = 1):
    
    #Check method functioning conditions
    if not isinstance(file, io.TextIOBase):
        raise ValueError('The object passed is not a text file object!')
    
    if not isinstance(num_line, int):
        raise ValueError('\'num_line\' parameter must be integer!')
    
    if num_col < 1:
        raise ValueError('\'num_col\' parameter must be greater than zero!')
    
    
    out_list = []
    n_line = 0
    
    #Read 'num_line' lines
    if num_line > 1: 
        
        for n in range(num_line):
            
            line = file.readline()
            n_line += 1
            
            #EOF error
            if line == '': 
                line = 'ERROR_read_list = EOF'
                return out_list, line, n_line
            
            tokens = line.split()
            #Inconsistent read error
            if len(tokens) < num_col:
                line = 'ERROR_read_list = Inconsistency'
                return out_list, line, n_line
            
            out_list.append(tokens)
        
        return out_list, line, n_line
            
            
    #Read lines consistently based on num_col
    elif num_line == 0:
        
        line = 'Start'
        while line != '':
            
            line = file.readline()
            n_line += 1
            #EOF may not be an error in this case
            if line == '': return out_list, line, n_line
            
            tokens = line.split()
            #Consistent reading condition
            if len(tokens) != num_col:
                break
            
            out_list.append(tokens)
        
        return out_list, line, n_line

###############################################################################

def errors_read_list(file_path, line, n_line):
    
    if 'EOF' in line:
        ErrorMSG = 'EOF reached while parsing file \'' + file_path \
            + '\' at line '+ str(n_line) +'!'
        raise ValueError(ErrorMSG)
        
    elif 'Inconsistency' in line:
        ErrorMSG = 'Inconsistent read while parsing the file \'' + file_path \
            + '\' at line '+ str(n_line) +'!'
        raise ValueError(ErrorMSG)
        
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
                    
                    
                    
                    