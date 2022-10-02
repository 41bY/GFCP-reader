# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 13:09:00 2022

@author: simula
"""
import numpy as np
# import os
import sys
from units import PhysConstants as phy
import Utils as ut

#Verify number of expected parameters
if len(sys.argv) != 4:
    if len(sys.argv) > 4:
        raise ValueError('Too many arguments')
        
    elif len(sys.argv) < 4:
        raise ValueError('Too few arguments')

#Get simulation input and output files and output directory, verify existance
IN_fname = str(sys.argv[1])
OUT_fname = str(sys.argv[2])
DATA_path = str(sys.argv[3])
ut.check_file_existance(IN_fname)
ut.check_file_existance(OUT_fname)
ut.create_dir(DATA_path)

#Create output files and open them to write
INFO_fname = DATA_path+'info.txt'
POS_fname = DATA_path+'Pos.xyz'
FUP_fname = DATA_path+'Forces_gf_up.dat'
FDW_fname = DATA_path+'Forces_gf_dw.dat'
TEM_fname = DATA_path+'Temp.dat'

f_POS = open(POS_fname, 'w')
f_FOR_up = open(FUP_fname, 'w')
f_FOR_up.writelines('Total forces on upper GF atoms: TimeStep(fs) Fx Fy Fz(GPa) \n')
f_FOR_dw = open(FDW_fname, 'w')
f_FOR_dw.writelines('Total forces on lower GF atoms: TimeStep(fs) Fx Fy Fz(GPa) \n')
f_TEM = open(TEM_fname, 'w')    
f_TEM.writelines('Temperature of GF atoms: TimeStep(fs) T(Kelvin) \n')  

#Read simulation input and initialize data for extraction
input_data = ut.read_CP_input(IN_fname)
iprint = input_data[0]
dt = input_data[1]
# celldm = input_data[2]
nat = input_data[3]
# ntyp = input_data[4]
cell = input_data[5]
types = input_data[6]

#Data extraction loop: read simulation output in blocks and write it to output files        
with open(OUT_fname, 'r') as f:
    step = -1
    expectation = None
    
    line = 'Start'
    while line != '':
            
        tags = ['kt_gf', 'upper[GPa]', 'lower[GPa]']
        tags_map, line = ut.read_tag(f, tags)
        if line == '': break
        
        line = ut.skip_lines(f, 1)
        if line == '': break
    
        if 'after' in line: #Collect data to write them, if 'after'
                    step += 1
                    line = ut.skip_lines(f, 1)
                    if line == '': break
                    
                    if 'ATOMIC_POSITIONS' in line: #Collect QM atom positions
                        qm_atoms, line = ut.read_list(f, num_line = nat, check_consistency = True)
                        if line == '': break
                        qm_atom_lbl = []
                        qm_atom_pos = []
                        for at in qm_atoms:
                            qm_atom_lbl.append(at[0])
                            qm_atom_pos.append(at[1:])
                        
                        qm_atom_lbl = np.array(qm_atom_lbl, dtype=str)
                        qm_atom_pos = np.array(qm_atom_pos, dtype=float)
                        
                    else:
                        expectation = '\'ATOMIC_POSITIONS\''
                        break    
                    
                    line = ut.skip_lines(f, nat+4)
                    if line == '': break
                    
                    if 'GF_ATOM_POSITIONS' in line: #Collect GF atom positions
                        gf_atoms, line = ut.read_list(f, num_line = 0, check_consistency = True)
                        if line == '': break
                        gf_atom_lbl = []
                        gf_atom_pos = []
                        for at in gf_atoms:
                            gf_atom_lbl.append('Cgf')
                            gf_atom_pos.append(at)
                        
                        gf_atom_lbl = np.array(gf_atom_lbl, dtype=str)
                        gf_atom_pos = np.array(gf_atom_pos, dtype=float)
                         
                    else:
                        expectation = '\'GF_ATOM_POSITIONS\''
                        break
                        
                    #Create all atom positions and labels
                    atom_pos = np.append(qm_atom_pos, gf_atom_pos, axis=0)
                    atom_lbl = np.append(qm_atom_lbl, gf_atom_lbl, axis=0)
                    
                    #At first iteration get slabs division and write '.info' file
                    if step == 0: 
                        nat_gf = gf_atom_pos[:,0].size
                        nat_tot = nat + nat_gf
                        mass_slabs = np.zeros(2, dtype=float) 
                    
                        atom_z = atom_pos[:,2]
                        idx_z_sort = np.argsort(atom_z)
                        atom_z = atom_z[idx_z_sort]
                        dz = atom_z[1:] - atom_z[:-1]
                        up_slab_index = np.argmax(dz) + 1
                        
                        #Calculate slab masses
                        for at in atom_lbl[idx_z_sort][:up_slab_index]:
                            if 'gf' in at: continue
                            mass_slabs[0] += types[at][1]
                        
                        for at in atom_lbl[idx_z_sort][up_slab_index:]:
                            if 'gf' in at: continue
                            mass_slabs[1] += types[at][1]
             
        
        elif 'refresh' in line: #Don't collect data if 'refresh'
                    line = ut.skip_lines(f, int(2*nat_tot + 7))
                    continue
        
        else:
                    expectation = '\'refresh\' or \'after\''
                    break
        
        #Skip irrelevant data
        line = ut.skip_lines(f, nat_gf)
        
        #Prepare data to write them to files
        time = float(step*dt)
        
        T_gf = tags_map['kt_gf']
        T_gf = T_gf.split()[0]
        T_gf = float(T_gf)
        
        F_gf_up = tags_map['upper[GPa]']
        F_gf_up = F_gf_up.split()
        F_gf_up = np.array(F_gf_up, dtype=float) 
        
        F_gf_dw = tags_map['lower[GPa]']
        F_gf_dw = F_gf_dw.split()
        F_gf_dw = np.array(F_gf_dw, dtype=float)
        
        #Order all atoms pos and lbl array w.r.t. z values at step 0
        atom_pos = atom_pos[idx_z_sort]
        atom_lbl = atom_lbl[idx_z_sort]
        
        #Pass from Bohr unit to angstrom
        atom_pos = atom_pos*phy.bohr_to_A
        
        #Write extracted data to files
        if step%1000 == 0: print('%d steps processed' %step)
        
        f_POS.writelines('%d \n' %nat_tot)
        f_POS.writelines('ATOMIC_POSITIONS (angstrom); Time = %15.10f \n' %time)
        for i, atom in enumerate(atom_pos):
            f_POS.writelines('%4s %25.15f %25.15f %25.15f \n' %(atom_lbl[i], atom[0], atom[1], atom[2]))
        
        f_FOR_up.writelines('%15.10f %25.15f %25.15f %25.15f \n' %(time, F_gf_up[0], F_gf_up[1], F_gf_up[2]))
        f_FOR_dw.writelines('%15.10f %25.15f %25.15f %25.15f \n' %(time, F_gf_dw[0], F_gf_dw[1], F_gf_dw[2]))
        f_TEM.writelines('%15.10f %25.15f \n' %(time, T_gf))
    
    types = list(types.values())
    step += 1 #Count the steps
    with open(INFO_fname, 'w') as f_info:
        
        f_info.writelines('Cell (angstrom):\n')
        for i in range(3):
            f_info.writelines('%15.10f %15.10f %15.10f' %(cell[i][0], cell[i][1], cell[i][2]))
            f_info.writelines('\n')
    
        f_info.writelines('\nAtom types (type|mass(a.m.u.)|#):\n')
        for at in types:
            f_info.writelines('%4s %10.4f %5d' %(at[0], at[1], at[2]))
            f_info.writelines('\n')
            
        f_info.writelines('\nTotal atoms number = %d \n' %nat_tot)
        f_info.writelines('\nQM atoms number = %d \n' %nat)
        f_info.writelines('\nGF atoms number = %d \n' %nat_gf)
        f_info.writelines('\nData time step = %10.9f fs\n' %(dt*iprint))
        f_info.writelines('\nSimulation time step = %10.9f fs\n' %(dt))
        
        f_info.writelines('\nData step = %d \n' %step)
        f_info.writelines('\nSimulation step (%dxData step) =  %d \n' %(step, step*int(iprint)))
        
        #NAT for down slab (the frist one) corresponds to starting index upper layer!
        f_info.writelines('\nSlabs info: (Slab \t #atoms \t Mass(a.m.u.))\n')
        f_info.writelines('Lower \t %d \t %10.5f \n' %(up_slab_index, mass_slabs[0]))
        f_info.writelines('Upper \t %d \t %10.5f \n' %(nat_tot - up_slab_index, mass_slabs[1]))    
    
    #Raise an error if the simulation output file is not as expected
    if expectation != None:
        f_POS.close()
        f_FOR_up.close()
        f_FOR_dw.close()
        f_TEM.close()
        raise ValueError('In simulation output file: \n line does not match the expectation \
                         %s in line: %s. \n File is corrupted or mismatch between simulation \
                              input and output %s %s!' %(expectation, line, IN_fname, OUT_fname))  

f_POS.close()
f_FOR_up.close()
f_FOR_dw.close()
f_TEM.close()  
            


# pos_old = []
# #Calculate all atom 's velocities
# if len(pos_old) == len(pos):
#     vel = []
#     for i, at in enumerate(pos):
#         Vx = ((at[1] - pos_old[i][1])/dt)*phy.AonFS_to_MonS
#         Vy = ((at[2] - pos_old[i][2])/dt)*phy.AonFS_to_MonS
#         Vz = ((at[3] - pos_old[i][3])/dt)*phy.AonFS_to_MonS
#         vel.append([at[0], Vx, Vy, Vz]) #m/s

# pos_old = pos
# write = True
            

# with open(INFO_fname, 'w') as f:
    
#     f.writelines('Cell (angstrom):\n')
#     for i in range(3):
#         f.writelines('%15.10f %15.10f %15.10f' %(cell[i][0], cell[i][1], cell[i][2]))
#         f.writelines('\n')

#     f.writelines('\nAtom types (type|mass(a.m.u.)|#):\n')
#     for at in types:
#         f.writelines('%4s %10.4f %5d' %(at[0], at[1], at[2]))
#         f.writelines('\n')
        
#     f.writelines('\nAtom number = %d \n' %nat)
#     f.writelines('\nData time step = %10.9f fs\n' %(dt*iprint))
#     f.writelines('\nSimulation time step = %10.9f fs\n' %(dt))               

# f_INFO = open(fname_INFO, 'a')
# f_INFO.writelines('\nData step: %d, (simulation step is %dxData step) \n' %(step, int(iprint)))
# f_INFO.writelines('\nSlabs info: (Slab \t Mass(a.m.u.))\n')
# f_INFO.writelines('Lower \t %10.5f \n' %M_slabs[0])
# f_INFO.writelines('Upper \t %10.5f \n' %M_slabs[1])

# f_INFO.writelines('\nOrdering POS, VEL for QM atoms: QM upper slab start at %d \n'\
#                   %i_split_QM)
# for i in i_ord_QM:
#     f_INFO.writelines('%4d ' %i)
    
# f_INFO.writelines('\nOrdering POS, VEL for GF atoms: GF upper slab start at %d \n'\
#                   %i_split_GF)
# for i in i_ord_GF:
#     f_INFO.writelines('%4d ' %i)

# f_INFO.close()            
                
                
                
                
                
                
                
                
            
            
            
