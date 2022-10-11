# -*- coding: utf-8 -*-
"""
Created on Sat Oct  1 10:22:10 2022

@author: simula
"""

###############################################################################
#Conversion factor

#Hartree time to femtoseconds
hatime2fs = 0.024189

#Bohr to angstrom
bohr2angstrom = 0.529177 

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
                    
                    
                    
                    