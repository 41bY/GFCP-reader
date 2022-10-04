# -*- coding: utf-8 -*-
"""
Created on Sat Oct  1 10:22:10 2022

@author: simula
"""
import numpy as np
import io
import UtilsLow as UtLow
import units

###############################################################################

def read_tag(file : io.TextIOBase, search_list : list, limit = np.infty):
    """
    Read 'file' based on tags contained in the list of strings 'search_list'.
    Once a tag is found in the file, its right hand side, delimited by '=', is read
    and stored within the output dictionary 'found_dict' whose keys correspond to
    the specified tags.

    Parameters
    ----------
    file : io.TextIOBase
        DESCRIPTION.
    search_list : list of str
        Tags to be searched and read in the file.
    limit : int, optional
        Number of line to read before stopping if some tag is not found.
        The default is np.infty i.e. read until EOF.

    Raises
    ------
    TypeError
        If 'file' is not a TextObject and if search_list is not a list of strings.
    ValueError
        If search_list is empty.

    Returns
    -------
    found_dict : dict
        The keys correspond to the specified tags, while the values are the strings
        read, associated to those tags. If some tags are not found before 'limit'
        or EOF, the values for those tags will be 'None'.
    line : str
        The last line read. If problems occur in the reading, 'line' will contain
        a string specifing the error: 
        'ERROR_read_tag = Limit reached', 'ERROR_read_tag = EOF', 'ERROR_read_tag = Corrupted field'
    n_line : int
        Number of lines read before return.

    """
    
    #Check method functioning conditions
    if not isinstance(file, io.TextIOBase):
        raise TypeError('The object passed is not a text file object!')

    if not isinstance(search_list, list):
        raise TypeError('The object passed is not a list!')        
    
    if len(search_list) == 0:
        raise ValueError('The search list cannot be empty!')
    
    for search in search_list:
        if not isinstance(search, str):
            raise TypeError('The search list must contain only strings!')
    
    found_dict = dict.fromkeys(search_list)
    
    line = 'Start'
    n_line = 0
    while line != '':
        
        #Limit the iterations
        if n_line > limit:
            line = 'ERROR_read_tag = Limit reached'
            return found_dict, line, n_line
        
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

def errors_read_tag(file_path : str, tags_map : dict, line : str, n_line : int):
    """
    Generate message errors refering to the method 'read_tag'. The generated 
    message errors depend on the value of 'line' in agreement with the method 'read_tag'.
    The parameter 'n_line' is used to specify at which line the reading error 
    is happened in the file given in 'file_path' while 'tags_map' is used to give
    further details on which tags have not be read.

    Parameters
    ----------
    file_path : str
        String containing the path of file processed by 'read_tag'
    tags_map : dict
        Dictionary containing the output dictionary returned by 'read_tag'.
    line : str
        String which specifies the type of errors encountered by 'read_tag'.
    n_line : int
        Number of the line at which the reading error is happened.

    Returns
    -------
    ErrorMSG : str
        Error message string.

    """
    
    if 'EOF' in line:
        not_found = [ele for ele in tags_map.keys() if tags_map[ele] == None]
        not_found_MSG = ''
        for nf in not_found: not_found_MSG += nf+' '
        
        ErrorMSG = 'EOF reached in file \'' + file_path + '\' at line '+ str(n_line) \
            + ': ' + not_found_MSG + 'not found!'
        return ErrorMSG
    
    elif 'Limit reached' in line:
        not_found = [ele for ele in tags_map.keys() if tags_map[ele] == None]
        not_found_MSG = ''
        for nf in not_found: not_found_MSG += nf+' '
        
        ErrorMSG = 'EOF reached in file \'' + file_path + '\' at line '+ str(n_line) \
            + ': ' + not_found_MSG + 'not found!'
        return ErrorMSG
        
    elif 'Corrupted field' in line:
        for tag, ele in tags_map.items():
            if ele == '':
                empty = tag
                break
            
        ErrorMSG = 'In file \'' + file_path + '\' ,corrupted field for tag: \'' + empty + \
            '\' at line ' + str(n_line)
        return ErrorMSG
    

###############################################################################

def read_list(file : io.TextIOBase, num_line : int = 0, num_col : int = 1):
    """
    Read 'num_line' lines sequentially in 'file' checking that they contain at 
    least 'num_col' columns after they have been splitted on blanks. If 'num_line' 
    is not given lines are read sequentually in a consistent manner i.e. each 
    read line must contain exactly the same number of columns specified by 'num_col'.
    
    Parameters
    ----------
    file : io.TextIOBase
        File to be read.
    num_line : int, optional
        Number of lines to be read. The default is 0.
    num_col : int, optional
        Minimum number of columns contained in the read lines. If 'num_line' is
        not given 'num_col' becomes the parameter for checking the consistency 
        of the read lines. The default is 1.

    Raises
    ------
    TypeError
        If 'file' is not a TextObject, if 'num_line' and 'num_col' are not integers.
    ValueError
        If 'num_col' is less than one i.e. at least 1 field must be specified
        and if 'num_line' is negative.

    Returns
    -------
    out_list : list of list of str
        The dimensions are 'num_line'x'num_col' if 'num_line' is specified while
        it is Nx'num_col' if it is not, where N is determined by the consistency.
    line : str
        The last line read. If problems occur in the reading, 'line' will contain
        a string specifing the error: 
        'ERROR_read_tag = EOF', 'ERROR_read_tag = Inconsistency'
    n_line : int
        Number of lines read before return.

    """
    #Check method functioning conditions
    if not isinstance(file, io.TextIOBase):
        raise TypeError('The object passed is not a text file object!')
    
    if not isinstance(num_line, int):
        raise TypeError('\'num_line\' parameter must be integer!')
    
    if not isinstance(num_col, int):
        raise TypeError('\'num_col\' parameter must be integer!')
    
    if num_line < 0:
        raise ValueError('\'num_line\' parameter cannot be negative!')
    
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

def errors_read_list(file_path : str, line : str, n_line : int):
    """
    Generate message errors refering to the method 'read_list'. The generated 
    message errors depend on the value of 'line' in agreement with the method 'read_list'.
    The parameter 'n_line' is used to specify at which line the reading error 
    is happened in the file given in 'file_path'.

    Parameters
    ----------
    file_path : str
        String containing the path of file processed by 'read_list'
    line : str
        String which specifies the type of errors encountered by 'read_list'.
    n_line : int
        Number of the line at which the reading error is happened.

    Returns
    -------
    ErrorMSG : str
        Error message string.

    """
    
    if 'EOF' in line:
        ErrorMSG = 'EOF reached while parsing file \'' + file_path \
            + '\' at line '+ str(n_line) +'!'
        return ErrorMSG
        
    elif 'Inconsistency' in line:
        ErrorMSG = 'Inconsistent read while parsing the file \'' + file_path \
            + '\' at line '+ str(n_line) +'!'
        return ErrorMSG
        
###############################################################################

def read_CP_input(file_path : str):
    """
    Method to read the input file of the 'Car-Parrinello' software (cp.x) in the 
    Quantum Espresso suite. The method is development to read only the part of 
    the file relevant for 'Extract_quantities.py'.
    It relies on the routines 'read_tag' and 'read_list' to read the file.

    Parameters
    ----------
    file_path : str
        Path to the 'Car-Parrinello' input file.

    Raises
    ------
    TypeError
        If 'file_path' is not a string.
    ValueError
        If even one tag specified in the returns of this method is not found in
        the CP input file, the method will raise a ValueError specifing which 
        tags are missing from the input.

    Returns
    -------
    iprint : int
        Number of dumped iteration in CP input file.
    dt : float
        Time step in CP input file.
    celldm : float
        Cell scale dimension in CP input file.
    nat : int
        Number of atoms in CP input file.
    ntyp : int 
        Number of atomic species in CP input file.
    cell : Numpy array of float
        Cell vectors in CP input file.
    types : dict
        Contain atomic information: keys are the different atomic species, 
        values are list containing (float) atomic mass and (int) number of atoms
        for each specie.

    """
    
    
    #Check method functioning conditions
    if not isinstance(file_path, str):
        raise TypeError('First argument must be a file path, it should be a str!')
    
    with open(file_path, 'r') as f:
        
        line = 'Start'
        n_line = 0
        namelist_found = 0
        while line != '':     
            
            line = f.readline()
            n_line += 1
            if line == '': 
                line = 'ERROR_eof'
                break
            
            if ('&CONTROL' in line) & (namelist_found == 0):
                tags = ['iprint', 'dt']
                tags_map, line, n_read = read_tag(f, tags)
                n_line += n_read
                if 'ERROR' in line: break
                
                iprint = tags_map['iprint']
                iprint = UtLow.trim_string(iprint, mode='num')
                #Check expected type
                if iprint == '':
                    line = 'ERROR = Corrupted field'
                    break
                iprint = int(iprint)
                
                dt = tags_map['dt']
                dt = UtLow.trim_string(dt, mode='num')
                #Check expected type
                if dt == '':
                    line = 'ERROR = Corrupted field'
                    break
                dt = float(dt)*2.4189*0.01 #fs
                
                namelist_found += 1
                    
            
            elif ('&SYSTEM' in line) & (namelist_found == 1):
                tags = ['celldm', 'nat', 'ntyp']
                tags_map, line, n_read = read_tag(f, tags)
                n_line += n_read
                if 'ERROR' in line: break
                
                celldm = tags_map['celldm']
                celldm = UtLow.trim_string(celldm, mode='num')
                #Check expected type
                if celldm == '':
                    line = 'ERROR = Corrupted field'
                    break
                celldm = float(celldm)
                

                nat = tags_map['nat']
                nat = UtLow.trim_string(nat, mode='num')
                #Check expected type
                if nat == '':
                    line = 'ERROR = Corrupted field'
                    break
                nat = int(nat)
                
                ntyp = tags_map['ntyp']
                ntyp = UtLow.trim_string(ntyp, mode='num')
                #Check expected type
                if ntyp == '':
                    line = 'ERROR = Corrupted field'
                    break
                ntyp = int(ntyp)
                
                namelist_found += 1
        
        
            elif ('CELL_PARAMETERS' in line) & (namelist_found == 2):
                token = line.split()
                
                #Read the unit length and convert it in Angstrom
                if len(token) > 1: #User defined unit length
                    unit = token[1]
                    
                    if 'bohr' in unit:
                        cell_unit = units.bohr_to_A
                    
                    elif 'alat' in unit:
                        cell_unit = celldm*units.bohr_to_A
                    
                    elif unit == 'angstrom':
                        cell_unit = 1.0
                        
                    else:
                        raise ValueError('Improper unit used in CELL_PARAMETERS\
                                         : %s' %unit)
                        
                else: #QE default unit length
                    cell_unit = celldm*units.bohr_to_A
                
                #Read the cell parameters              
                cell, line, n_read = read_list(f, num_line = 3, num_col = 3)
                n_line += n_read
                if 'ERROR' in line: break
                cell = cell_unit*np.array(cell, dtype=float)
                
                namelist_found += 1
                
                
            elif ('ATOMIC_SPECIES' in line) & (namelist_found == 3):
                atom_species, line, n_read = read_list(f, num_line = ntyp, num_col = 3)
                n_line += n_read
                if 'ERROR' in line: break
                
                #Create atom types list: [kind, mass, number]
                atom_kinds = [(at[0], [at[0], float(at[1]), 0]) for at in atom_species]
                types = dict(atom_kinds)
                
                namelist_found += 1
    
                
            elif ('ATOMIC_POSITIONS' in line) & (namelist_found == 4):
                atom_pos, line, n_read = read_list(f, num_line = nat, num_col = 4)
                n_line += n_read
                if 'ERROR' in line: break
                
                #Counts the number of each different atom and insert its number
                for at in atom_pos:
                    types[at[0]][2] += 1                  
                
                #Last readlist, when completed terminate the reading
                break
        
        if 'ERROR' in line:
            
            if '_read_tag' in line:
                ErrorMSG = errors_read_tag(file_path, tags_map, line, n_line)
                raise ValueError(ErrorMSG)
                
            elif '_read_list' in line:
                ErrorMSG = errors_read_list(file_path, line, n_line)
                raise ValueError(ErrorMSG)
                
            elif '_eof' in line:
                raise ValueError('EOF reached at line %d: only %d namelists found out of the specified 5 namelists:\
                                 &CONTROL &SYSTEM CELL_PARAMETERS ATOMIC_SPECIES ATOMIC_POSITIONS'\
                                     %(n_line, namelist_found))
            
            elif 'Corrupted field' in line:
                ErrorMSG = 'In file \'' + file_path + '\' ,corrupted field at line ' + str(n_line)
                raise ValueError(ErrorMSG)
            
            else:
                raise ValueError('Unaxpected error!')
                

    return iprint, dt, celldm, nat, ntyp, cell, types                           
                    
###############################################################################

def read_ref_output(file_path : str):
    """
    Reference read of the output file for the GFCP software (gfcp.x).

    Parameters
    ----------
    file_path : str
        File name path.

    Raises
    ------
    TypeError
        If 'file_path' is not a string.

    Returns
    -------
    atom_label : list of str
        Containing the atomic specie labels of all the atoms.
    atom_pos : Nx3 Numpy array of float
        Containing the reference spatial coordinates of all the atoms.
    line : str
        Last read line of the file.
    n_line : int
        Number of lines read.

    """
    #Check method functioning conditions
    if not isinstance(file_path, str):
        raise TypeError('First argument must be a file path, it should be a str!')
    
    with open(file_path, 'r') as f:
        
        line = 'Start'
        n_line = 0
        
        #Save qm atoms and gf atoms as soon as possible
        while line != '':
            
            line, n_read = UtLow.skip_lines(f, 1)
            n_line += n_read
            if 'ATOMIC_POSITIONS' in line:
                #Collect QM atom positions
                qm_atoms, line, n_read = read_list(f, num_line = 0, num_col = 4)
                n_line += n_read
                if 'ERROR' in line: return line, n_line
                
                while line != '':
                    
                    line, n_read = UtLow.skip_lines(f, 1)
                    n_line += n_read
                    if 'GF_ATOM_POSITIONS' in line: 
                        #Collect GF atom positions
                        gf_atoms, line, n_read = read_list(f, num_line = 0, num_col = 3)
                        n_line += n_read
                        if 'ERROR' in line: return line, n_line
                        
                        break
                    
                    elif 'ERROR' in line: return line, n_line
                break
            
            elif 'ERROR' in line: return line, n_line
            
    qm_atom_lbl = []
    qm_atom_pos = []
    for at in qm_atoms:
        qm_atom_lbl.append(at[0])
        qm_atom_pos.append(at[1:])
    
    qm_atom_pos = np.array(qm_atom_pos, dtype=float)
    
    gf_atom_lbl = []
    gf_atom_pos = []
    for at in gf_atoms:
        gf_atom_lbl.append('Cgf')
        gf_atom_pos.append(at)
    
    gf_atom_pos = np.array(gf_atom_pos, dtype=float)
    
    atom_label = qm_atom_lbl + gf_atom_lbl
    atom_pos = np.append(qm_atom_pos, gf_atom_pos, axis=0)
    
    return atom_label, atom_pos, line, n_line

###############################################################################        
    
def compatibility_check(labels_ref : list, atom_types : dict):
    """
    Check that the CP input and GFCP output are compatible i.e. that the atomic
    systems defined in CP input and GFCP output are the same. This is acomplished
    by comparing the two dictionaries specifing the atomic systems.

    Parameters
    ----------
    labels_ref : list of str
        Containing the atomic specie labels of all the atoms found in the output
        file by the method 'read_ref_output'. This list is transformed into a
        dictionary by this method.
    atom_types : dict
        Contain atomic information returned from the reading of CP input by the
        method 'read_CP_input': keys are the different atomic species, 
        values are list containing (float) atomic mass and (int) number of atoms
        for each specie. 

    Raises
    ------
    TypeError
        If 'labels_ref' is not a list of str and if 'atom_types' is not a dict
        of list as returned from the method 'read_CP_input'.

    Returns
    -------
    None.

    """
    
    #Check method functioning conditions
    if not isinstance(labels_ref, list):
        raise TypeError('\' labels_ref \' is not a list!')
    
    if not isinstance(atom_types, dict):
        raise TypeError('\' atom_types \'is not a dictionary!')
    
    #Build types_ref dictionary
    labels_ref = np.array(labels_ref, dtype=str)
    lbl_ref, num_ref = np.unique(labels_ref, return_counts=True)
    types_ref = dict([[lbl, num_ref[i]] for i, lbl in enumerate(lbl_ref)])
    
    #Build types dictionary without the atom masses
    types = {}
    for key in atom_types.keys():
        types[key] = atom_types[key][2]
        
    if types_ref == types: return True
    else: return False

###############################################################################

def get_slab_info(atom_lbl, atom_pos, atom_types):
    """
    Get slabs information from returns of the reference read of GFCP output 
    performed by the method 'read_ref_output'.

    Parameters
    ----------
    atom_lbl : list of str
        Containing the atomic specie labels of all the atoms.
    atom_pos : Nx3 Numpy array of float
        Containing the reference spatial coordinates of all the atoms.
    atom_types : dict
        Contain atomic information returned from the reading of CP input by the
        method 'read_CP_input': keys are the different atomic species, 
        values are list containing (float) atomic mass and (int) number of atoms
        for each specie. 

    Returns
    -------
    mass_slabs : list of float
        Contain the masses of each slabs.
    up_slab_index : int
        Index of the first upper slab atom in the sorted 'atom_lbl' and 'atom_pos'.
    idx_z_sort : list of int
        Sorted indeces for 'atom_lbl' and 'atom_pos' by their z coordinate 
        ascending value.

    """

    #Get slab separation---------------------------------------------------
    atom_z = atom_pos[:,2]
    idx_z_sort = np.argsort(atom_z)
    atom_z = atom_z[idx_z_sort]
    dz = atom_z[1:] - atom_z[:-1]
    up_slab_index = np.argmax(dz) + 1
    
    mass_slabs = np.zeros(2, dtype=float) 
    #Calculate slab masses-------------------------------------------------
    atom_lbl = np.array(atom_lbl, dtype=str)
    for at in atom_lbl[idx_z_sort][:up_slab_index]:
        if 'gf' in at: continue
        mass_slabs[0] += atom_types[at][1]
    
    for at in atom_lbl[idx_z_sort][up_slab_index:]:
        if 'gf' in at: continue
        mass_slabs[1] += atom_types[at][1]
    
    return mass_slabs, up_slab_index, idx_z_sort
    
###############################################################################

def write_xyz(file : io.TextIOBase, data : tuple, header : str):
    """
    Write an '.xyz' formatted file. The format is:
        N
        Comment line
        atom_1_label     x1  y1  z1
        atom_2_label     x2  y2  z2
        ...                 ...
        atom_N_label     xN  yN  zN
        
    Parameters
    ----------
    file : io.TextIOBase
        File to be written.
    data : tuple
        Contains a N dimensional list of str representing the atomic label 
        as first element and a Nx3 Numpy array representing the atomic coordinates
        as second element.
    header : str
        Comment line present in the '.xyz' format.

    Raises
    ------
    TypeError
        If 'file' is not a text file object and if 'data' is not a tuple
    ValueError
        If 'data[0]' and 'data[1]' do not have the same dimension along the 
        first axis.

    Returns
    -------
    None.

    """
    #Check method functioning conditions
    if not isinstance(file, io.TextIOBase):
        raise TypeError('\'file\' is not a text file object!')
    
    if not isinstance(data, tuple):
        raise TypeError('\'data\' is not a tuple!')
    
    if not isinstance(header, str):
        raise ValueError('\'header\' is not a string!')
    
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
        raise TypeError('\'file_path\' must be a file path, it should be a str!')
    
    if not isinstance(info_data, list):
        raise TypeError('Info_data must be a list!')
    
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
