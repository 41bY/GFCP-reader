# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 13:09:00 2022

@author: simula
"""
import numpy as np
import os
import sys
import Read_input as ri
import Read_reference_output as rr
import Read_output as ro
import Utils as ut

#Command line managment section------------------------------------------------
#Verify number of expected parameters
if len(sys.argv) != 4:
    if len(sys.argv) > 4:
        raise ValueError('Too many arguments: 3 expected!')
        
    elif len(sys.argv) < 4:
        raise ValueError('Too few arguments: 3 expected!')

#Get simulation input, simulation output and output directory
IN_fname = str(sys.argv[1])
OUT_fname = str(sys.argv[2])
DATA_path = str(sys.argv[3])

#Verify existance of simulation files and create output directory if not present
if not os.path.isfile(IN_fname): raise ValueError('%s not found' %IN_fname)
if not os.path.isfile(IN_fname): raise ValueError('%s not found' %OUT_fname)
if not os.path.isdir(DATA_path): os.mkdir(DATA_path)

#Extraction initialization section---------------------------------------------
#Read simulation input and initialize data for extraction
iprint, dt, celldm, nat_qm, ntyp, cell, types, label_qm_in = ri.read_CP_input(IN_fname)

#Reference reading of simulation output
tags_line_num, block_size, label_qm_out, nat_gf, idx_z_sort, idx_up_slab = rr.read_ref_output(OUT_fname)

#Compatibility check between simulation input and output
if label_qm_in != label_qm_out: raise ValueError('Simulation input and output are not compatible!')

#Output files initialization
f_TEM, f_FOR_up, f_FOR_dw, f_POS = ro.initialize_output_files(DATA_path)

#Data extraction section-------------------------------------------------------        
with open(OUT_fname, 'r') as file:
    
    #Skip comment lines
    for i in range(3): next(file)
    
    #Initialize error, number of relax and iteration
    err = ''
    num_relax = 0
    iteration = 0
    
    #Iteration loop: read simulation output in iteration blocks, check for errors and write data to output files
    while err == '':
        
        #Read iteration block from data file
        block_list = ro.read_block(file, block_size)
        
        #Exit the loop if EOF
        if block_list[0] == '': break
        
        #Slice block to get temperature, forces, control, qm atoms, gf atoms
        temp_str, fup_str, fdw_str, write_str, qm_list, gf_list = ro.get_str_values(block_list, tags_line_num, nat_qm, nat_gf)
        
        #write_str determines if data has to be written or not
        if 'after' in write_str: pass
    
        elif 'relax' in write_str: 
            num_relax += 1
            continue
        
        else: #Bad reading
            err = 'Expected \'after\' or \'relax\' at line'
            err += ro.Error_msg(iteration, num_relax, block_size, tags_line_num[3])
            break
        
        #Check and organize data, if corrupted fields throw errors
        #Temperature data
        temp, err = ro.check_temp(temp_str)
        if err != '': 
            err += ro.Error_msg(iteration, num_relax, block_size, tags_line_num[0])
            break
        
        #Force up data
        fup, err = ro.check_force(fup_str)
        if err != '': 
            err += ro.Error_msg(iteration, num_relax, block_size, tags_line_num[1])
            break
        
        #Force down data
        fdw, err = ro.check_force(fdw_str)
        if err != '':
            err += ro.Error_msg(iteration, num_relax, block_size, tags_line_num[2])
            break
        
        #QM atoms data
        qm_atoms, err, n_line = ro.check_qm_atoms(qm_list)
        if err != '':
            err += ro.Error_msg(iteration, num_relax, block_size, n_line+tags_line_num[4])
            break
        
        #GF atoms data
        gf_atoms, err, n_line = ro.check_gf_atoms(gf_list)
        if err != '': 
            err += ro.Error_msg(iteration, num_relax, block_size, n_line+tags_line_num[5])
            break     
    
        #Order all atoms as in reference
        atoms = ro.order_all_atoms(qm_atoms, gf_atoms, idx_z_sort)
        
        #Prepare data format to write them to output files
        time = iteration*dt*iprint*ut.hatime2fs
        tempVStime, fupVStime, fdwVStime, atom_label, atom_pos, header = ro.get_formatted_data(time, temp, fup, fdw, atoms)
        atom_pos = atom_pos*ut.bohr2angstrom
        
        #Write output files
        np.savetxt(f_TEM, tempVStime, fmt=['%10.5f', '%25.12f'])
        np.savetxt(f_FOR_up, fupVStime, fmt=['%10.5f', '%25.12f', '%25.12f', '%25.12f'])
        np.savetxt(f_FOR_dw, fupVStime, fmt=['%10.5f', '%25.12f', '%25.12f', '%25.12f'])
        ro.write_xyz(f_POS, atom_label, atom_pos, header)
        
        #Update iteration counter
        iteration += 1
    
    #If iteration loop stops, close output files
    f_POS.close()
    f_FOR_up.close()
    f_FOR_dw.close()
    f_TEM.close()    
    
    #Handle errors
    if err != '':
        err_msg = err + ' Successfully read ' + str(iteration) + ' iterations'
        raise ValueError(err_msg)
    
    #No-error read
    else:
        msg = 'Read of simulation output completed. ' + str(iteration) +' iterations read successfully.'
        print(msg)