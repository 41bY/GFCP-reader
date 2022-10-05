# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 11:10:31 2022

@author: simula
"""
import UtilsLow as UtLow
import UtilsHigh as UtHig
import os


#-----------------------------------------------------------------------------#
#---------------------------UtilsLow Tests------------------------------------#
#-----------------------------------------------------------------------------#
###############################################################################

def test_trim_string():
    numbers = '123'
    alphas = 'abc'
    
    for c in numbers:
        expected = c in UtLow.trim_string('abc123', mode='alpha')
        assert (expected == False)
        
    for c in alphas:
        expected = c in UtLow.trim_string('abc123', mode='alpha')
        assert (expected == True)
        
###############################################################################

def test_check_file_existance():
    exist = './'
    not_exist = './notexist.txt'
    
    assert UtLow.check_file_existance(exist)
    assert UtLow.check_file_existance(not_exist)

###############################################################################

def test_create_dir():
    TESTDIR = './Testing/TESTDIR'
    
    assert TESTDIR == UtLow.create_dir(TESTDIR)
    assert os.path.isdir(TESTDIR)

###############################################################################
    
def test_skip_lines():
    file = open('./Testing/TESTDIR/File_to_test_routines.txt', 'r')
    assert UtLow.skip_lines(file, 4)[0] == 'Line 4: \n'
    file.close()
    
    file = open('./Testing/TESTDIR/File_to_test_routines.txt', 'r')
    assert UtLow.skip_lines(file, 4)[1] == 4
    file.close()
    
    file.close()
    
    file = open('./Testing/TESTDIR/File_to_test_routines.txt', 'r')
    assert UtLow.skip_lines(file, 10)[0] == 'ERROR = EOF'
    file.close()
    
    file = open('./Testing/TESTDIR/File_to_test_routines.txt', 'r')
    assert UtLow.skip_lines(file, 10)[1] == 7
    file.close()
    
###############################################################################
    
    
    
    
    
    
    
    
    
    
    