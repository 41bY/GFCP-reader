# -*- coding: utf-8 -*-
"""
Created on Sat Oct  1 10:22:10 2022

@author: simula
"""

import numpy as np
import Utils as ut

###############################################################################

def read_CP_input(file_path : str):
    """
    'main' program of this module. It reads a 'cp.x' input file and it controls
    the other methods above. Open the 'cp.x' input file, load it into list of str,
    call 'search_input_tags' to get the string values associated to the tags and
    perform checks on them by calling 'check_founds'.
    It returns the list of founds associated to the tags.

    Parameters
    ----------
    file_path : str
        'cp.x' input file path.

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

    """

    #Input files are small so they can be load into memory
    file = open(file_path, 'r')
    lines = file.readlines()
    file.close()
    
    #Search tags in lines and returns string values inside found_list
    found_list_str = search_input_tags(lines)
    
    #Check formats -> raise errors if tag not found or corrupted. Assign right types to list values
    found_list = check_founds(found_list_str)
    
    return found_list

###############################################################################

def search_input_tags(read_list : list):
    """
    Search tags needed to perform the data extracion. The searching is performed
    within 'read_list' containing the content, line per line, of the 'cp.x' input
    file. It returns a list containing the string associated to each tag:
    'iprint', 'dt', 'celldm', 'nat', 'ntyp', 'CELL_PARAMETERS', 'ATOMIC_SPECIES',
    'ATOMIC_POSITIONS'

    Parameters
    ----------
    read_list : list
        Input list containing the lines of 'cp.x' input file.

    Returns
    -------
    found_list : list of str
        Output list containing the str values associated to each tag. 
        If a tag is corrupted the associated str will be ''. 
        If a tag is missing 'found_list' will have a smaller size.
    """
    
    #Initialize variables expected to be found in CP input file
    found_list = []
    nat = 0
    ntyp = 0
    
    #Search loop
    for num_line, line in enumerate(read_list):
        
        if 'iprint' in line:
            value_str = process_tag(line)
            found_list.append(value_str)
            
        elif 'dt' in line:
            value_str = process_tag(line)
            found_list.append(value_str)
            
        elif 'celldm' in line:
            value_str = process_tag(line)
            found_list.append(value_str)
            
        elif 'nat' in line:
            value_str = process_tag(line)
            found_list.append(value_str)
            if value_str != '': nat = int(value_str)
            
        elif 'ntyp' in line:
            value_str = process_tag(line)
            found_list.append(value_str)
            if value_str != '': ntyp = int(value_str)
            
        elif 'CELL_PARAMETERS' in line:
            cell = process_cell(read_list[num_line:])
            found_list.append(cell)

        elif 'ATOMIC_SPECIES' in line:
            species = process_species(read_list[num_line:], ntyp)
            found_list.append(species)

        elif 'ATOMIC_POSITIONS' in line:
            pos = process_atoms(read_list[num_line:], nat)
            found_list.append(pos)                      
            
    return found_list

###############################################################################

def process_tag(input_str : str):
    """
    Standard way to process a line associated to a tag. The line is splitted
    with '=' and the right-hand side is trimmed to contain only numbers.
    This processed string is returned. If the field is corrupted a null string
    is returned.

    Parameters
    ----------
    input_str : str
        Input string associated to a tag.

    Returns
    -------
    value_str : str
        The trimmed righ-hand side of the line. 
        Null string if field is corrupted.

    """
    
    tokens = input_str.split('=')
    value_str = ''
    
    #Error, return blanck string
    if len(tokens) < 2:
        value_str
        return value_str
        
    else:
        token = tokens[1]
        value_str = ut.trim_string(token, mode = 'num')
        
        #Error, return blanck string
        if value_str == '':
            return value_str
        
        else:
            return value_str

###############################################################################

def process_cell(input_list_str : list):
    """
    Standard way to process a 'CELL_PARAMTERS' tag in 'cp.x' input file. 
    Three lines below the tag are read, splitted and trimmed to create a 3x3
    matrix of string containing the lattice vectors.

    Parameters
    ----------
    input_list_str : list
        Input list of string with length 4. Containing the tag and the three lattice
        vectors.

    Returns
    -------
    cell : list of str
        If the field is not corrupted it returns a 3x3 list of str containing the
        cartesian component of the three lattice parameters.
        If the field is corrupted it returns a void list.

    """
    
    cell = []
    
    #Error return void cell
    if len(input_list_str) < 4:
        return cell
    
    else:
        
        for line in input_list_str[1:4]:
            tokens = line.split()
            
            #Error return void cell
            if len(tokens) != 3:
                cell = []
                return cell
            
            else:
                row = []
                for tk in tokens:
                    token = ut.trim_string(tk, mode = 'num')
                    
                    #Error return void cell
                    if token == '':
                        cell = []
                        return cell
                    
                    else:
                        row.append(token)

                cell.append(row)
                
        return cell

###############################################################################

def process_species(input_list_str : list, ntyp : int):
    """
    Standard way to process a 'ATOMIC_SPECIES' tag in 'cp.x' input file. 
    'ntyp' lines below the tag are read, splitted and trimmed to create a 'ntyp'x2
    matrix of string containing the different atomic labels and atomic masses.
    'ntyp' has to be provided.

    Parameters
    ----------
    input_list_str : list
        Input list of str of dimension ('ntyp'+1)x3. The first row contains the
        tag 'ATOMIC_SPECIES'. The other rows contain atomic labels atomic masses 
        and atomic pseudopotential.
    ntyp : int
        Number of different atomic species present in the 'cp.x' input file.

    Returns
    -------
    species : list
        List of str of dimension 'ntyp'x2 containing string of atomic labels 
        and masses. If the tag is corrupted a void list is returned.

    """
    
    
    species = []
    
    #Error, returns void species
    if len(input_list_str) < ntyp + 1:
        return species
    
    else:
        
        for line in input_list_str[1:ntyp + 1]:
            tokens = line.split()
            
            #Error, returns void species
            if len(tokens) < 3:
                species = []
                return species
            else:
                typ = tokens[0]
                token = tokens[1]
                mass_str = ut.trim_string(token, mode = 'num')
                
                #Error, return void species
                if mass_str == '':
                    species = []
                    return species
                    
                else:
                    row = [typ, mass_str]

                species.append(row)

        return species

###############################################################################

def process_atoms(input_list_str : list, nat : int):
    """
    Standard way to process a 'ATOMIC_POSITIONS' tag in 'cp.x' input file. 
    'nat' lines below the tag are read, splitted and trimmed to create a 'nat'x1
    matrix of string containing the atomic labels of each atom in the input.
    'nat' has to be provided.

    Parameters
    ----------
    input_list_str : list
        Input list of str of dimension ('nat'+1)x4. The first row contains the
        tag 'ATOMIC_POSITIONS'. The other rows contain atomic labels and the 
        three cartesian component of each atom in input.
    nat : int
        Total number of atoms in the 'cp.x' input file.

    Returns
    -------
    atoms : list
        List of str of dimension 'nat' containing string of atomic labels for
        each atom. If the tag is corrupted a void list is returned.

    """
    
    
    atoms = []
    
    #Error return void list
    if len(input_list_str) < nat + 1:
        atoms = []
        return atoms
    
    else:
        
        for line in input_list_str[1:nat + 1]:
            tokens = line.split()
            
            #Error return void list
            if len(tokens) < 1:
                atoms = []
                return atoms
            
            else:
                typ = tokens[0]
                atoms.append(typ)

        return atoms

###############################################################################

def check_founds(found_list : list):
    """
    Check 'found_list' has the proper length i.e. every tag has been found and
    verify that each tag is not void or blanck. In this case they are converted
    to their proper type. Otherwise an exception occurs specifing that the 
    'cp.x' input file is corrupted.

    Parameters
    ----------
    found_list : list of str
        Input list containing the str values associated to each tag. 
        If a tag is corrupted the associated str will be ''. 
        If a tag is missing 'found_list' will have a smaller size.

    Raises
    ------
    ValueError
        Tag not found if 'found_list' has a different dimension that the one expected.
        Corrupted tag if a specific tag is blanck.

    Returns
    -------
    found_list : list
        Output list containing the values associated to each tag of the proper
        type depending on the tag: 
        'iprint' -> int; 'dt' -> float; 'celldm' -> float; 'nat' -> int; 'ntyp' -> int;
        'CELL_PARAMETERS' -> np.array(float); 'ATOMIC_SPECIES' -> list of [str, float];
        'ATOMIC_POSITIONS' -> list of str

    """
    
    if len(found_list) != 8:
        raise ValueError('Tags: \'iprint\' | \'dt\' | \'celldm\' | \'ntyp\' |'\
                         +' \'nat\' | \'CELL_PARAMETERS\' | \'ATOMIC_POSITIONS\' | \'ATOMIC_SPECIES\''\
                         +' not found in CP input file.')
    
    #iprint
    if found_list[0] != '':
        found_list[0] = int(found_list[0])
    else:
        raise ValueError('Corrupted tag \'iprint\' in CP input file')
    
    #dt
    if found_list[1] != '':
        found_list[1] = float(found_list[1])
    else:
        raise ValueError('Corrupted tag \'dt\' in CP input file')
    
    #celldm
    if found_list[2] != '':
        found_list[2] = float(found_list[2])
    else:
        raise ValueError('Corrupted tag \'celldm\' in CP input file')

    #nat
    if found_list[3] != '':
        found_list[3] = int(found_list[3])
    else:
        raise ValueError('Corrupted tag \'nat\' in CP input file')

    #ntyp
    if found_list[4] != '':
        found_list[4] = int(found_list[4])
    else:
        raise ValueError('Corrupted tag \'ntyp\' in CP input file')
        
    #Cell
    if found_list[5] != []:
        found_list[5] = np.array(found_list[5], dtype=float)
    else:
        raise ValueError('Corrupted tag \'CELL_PARAMETERS\' in CP input file')
    
    #Species
    if found_list[6] != []:
        found_list[6] = [[fnd[0], float(fnd[1])] for fnd in found_list[6]]
    else:
        raise ValueError('Corrupted tag \'ATOMIC_SPECIES\' in CP input file')
    
    #Species
    if found_list[7] != []:
        found_list[7] = found_list[7]
    else:
        raise ValueError('Corrupted tag \'ATOMIC_POSITIONS\' in CP input file')
        
    return found_list

###############################################################################