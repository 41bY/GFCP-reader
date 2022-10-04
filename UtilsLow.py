# -*- coding: utf-8 -*-
"""
Created on Sat Oct  1 10:22:10 2022

@author: simula
"""

import os
import io

###############################################################################

def trim_string(in_str : str, mode = 'alpha-num'):
    """
    Trim a string depending on 'mode' values. Default 'mode' is set to 'alpha-num'
    which stands for alpha-numeric. Other default implementation is 'alpha' which
    stands for alphabetic and 'num' numeric. Eventually one can give a user defined
    mode by specifing what to include in the trimming process.

    Parameters
    ----------
    in_str : str
        Input string to be trimmed.
    mode : str, optional
        Trimming mode. The default is 'alpha-num' and stands for numeric trimming
        where letters, numbers and '.' are conserved in the string. The other default
        implementation is 'alpha' which stands for alphabetic where only letters
        are conserved in the string. By passing a generic string to mode the
        trimming will be personalized based on that string.

    Returns
    -------
    out_str : str
        The trimmed string derived from 'in_str' based on 'mode'.

    """
    if mode == 'alpha-num':
        out_str = ''.join([c for c in in_str if c in '0123456789.abcdefghkjilmnopqrstuvwxyz'])    
    
    elif mode == 'num':
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
    """
    Skip 'lines_to_skip' number of lines in 'file'. The pointer to the line file
    is parsed 'lines_to_skip' lines ahead of its initial position.

    Parameters
    ----------
    file : io.TextIOBase
        File where the lines have to be skipped. 
    lines_to_skip : int
        Number of lines to be skipped.

    Raises
    ------
    TypeError
        Raise a TypeError if 'file' is not a text object file and if 'lines_to_skip'
        is not an integer.
    
    ValueError
        Raise a ValueError if 'lines_to_skip' is not a positive integer.
    Returns
    -------
    line : str
        The last skipped line. If the skipped line is EOF, 'line' will be a string
        containing 'ERROR = EOF'
    n_line : int
        Number of lines skipped before return.

    """
    
    #Check method functioning conditions
    if not isinstance(file, io.TextIOBase):
        raise TypeError('The object passed is not a text file object!')
        
    if not isinstance(lines_to_skip, int):
        raise TypeError('Number of lines to be skipped must be an integer!')        
    
    if not lines_to_skip > 0:
        raise ValueError('Number of lines to be skipped must be positive!')
    
    n_line = 0
    for i in range(lines_to_skip):
        line = file.readline()
        n_line += 1
        
        if line == '':
            line = 'ERROR = EOF'
            return line, n_line
    
    return line, n_line

###############################################################################