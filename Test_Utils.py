# -*- coding: utf-8 -*-
"""
Created on Sat Oct  1 10:22:10 2022

@author: simula
"""

import Utils as ut

###############################################################################

def test_trim_string_alphanum():
    """
    Test the success of the function to trim a string on the default mode 'alpha-num'.

    Returns
    -------
    None.

    """
    
    sample = 'ab!cd3??456erd21#@'

    trimmed = ut.trim_string(sample, mode = 'alpha-num', trim_from = '')
    
    assert trimmed == 'abcd3456erd21'

###############################################################################

def test_trim_string_alpha():
    """
    Test the success of the function to trim a string on the default mode 'alpha'.

    Returns
    -------
    None.

    """
    
    sample = 'ab!cd3??456erd21#@'

    trimmed = ut.trim_string(sample, mode = 'alpha', trim_from = '')
    
    assert trimmed == 'abcderd'

###############################################################################

def test_trim_string_num():
    """
    Test the success of the function to trim a string on the default mode 'num'.

    Returns
    -------
    None.

    """
    
    sample = 'ab!cd3??456erd21#@'

    trimmed = ut.trim_string(sample, mode = 'num', trim_from = '')
    
    assert trimmed == '345621'

###############################################################################

def test_trim_string_user():
    """
    Test the success of the function to trim a string on the user defined 'trim_from'

    Returns
    -------
    None.

    """
    
    sample = 'ab!cd3??456erd21#@'

    trimmed = ut.trim_string(sample, mode = 'user', trim_from = '#@?!')
    
    assert trimmed == '!??#@'

###############################################################################        
                    
                    
                    
                    