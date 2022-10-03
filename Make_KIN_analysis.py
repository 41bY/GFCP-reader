# -*- coding: utf-8 -*-
"""
Created on Sat Sep 10 16:09:07 2022

@author: simula
"""
import numpy as np
import sys
import Utils as ut
from units import PhysConstants as phy

#Verify number of expected parameters
if len(sys.argv) != 4:
    if len(sys.argv) > 4:
        raise ValueError('Too many arguments')
        
    elif len(sys.argv) < 4:
        raise ValueError('Too few arguments')
        
#Get simulation input and output files and output directory, verify existance
INFO_fname = str(sys.argv[1])
POS_fname = str(sys.argv[2])
DATA_path = str(sys.argv[3])
ut.check_file_existance(INFO_fname)
ut.check_file_existance(POS_fname)
ut.create_dir(DATA_path)

#Create output files
VEL_fname = DATA_path+'Vel.xyz'
VCOM_fname = DATA_path+'VCOM.dat'
RCOM_fname = DATA_path+'RCOM.dat'

#Open output files to write them in the main loop
f_VEL = open(VEL_fname, 'w')
f_VCOM = open(VCOM_fname, 'w')
f_VCOM.writelines('Relative slabs center of mass velocity: TimeStep(fs) Vx Vy Vz(m/s) \n')
f_RCOM = open(RCOM_fname, 'w')
f_RCOM.writelines('Relative slabs center of mass position: TimeStep(fs) Rx Ry Rz(angstrom) \n')

#Load atomic info
with open(INFO_fname, 'r') as f:

    line = 'Start'
    while line != '':
        
        line = f.readline()
        
        if 'Atom types' in line:
            atoms, line = ut.read_list(f, num_line = 0, check_consistency = True)
            if line == '': break
            
            types = dict([[atom[0], float(atom[1])] for atom in atoms])

        elif 'Slabs info' in line:
            slabs, line = ut.read_list(f, num_line = 2, check_consistency = True)
            if line == '': break
            
            up_slab_index = int(slabs[0][1])
            mass_slabs = [float(slab[2]) for slab in slabs]

            
            
#Create 'Vel.xyz', 'COM.dat' and 'VCOM.dat' using 'Pos.xyz'
with open(POS_fname, 'r') as f:
    
    #Initialize main variables
    line = 'Start'
    step = -1
    while line != '':
        
        line = f.readline()
        if line == '': break
        nat = int(line) #Number of atoms is expected in a '.xyz'
        
        line = f.readline()
        if line == '': break
        time = float(line.split('=')[1]) #Number of atoms is expected in a '.xyz'        

        #Read atom positions and create label and pos arrays
        atoms, line = ut.read_list(f, num_line = nat, check_consistency = True)
        at_lbl = []
        at_pos = []
        for atom in atoms:
            at_lbl.append(atom[0])
            at_pos.append(atom[1:])
        at_pos = np.array(at_pos, dtype=float)   
            
        #At the first iteration: calculate weigths array to obtain RCOM and VCOM
        #Save position and times to calculate velocity at the next next iteration
        if step < 0:
            at_pos_2old = at_pos
            time_2old = time
            
            n_gf_slabs = [0, 0]
            weigths = []
            for i, lbl in enumerate(at_lbl):
                if i < up_slab_index:
                    if 'gf' in lbl:
                        n_gf_slabs[0] += 1
                        continue
                    else:
                        weigth = -types[lbl]/mass_slabs[0]
                
                else:
                    if 'gf' in lbl:
                        n_gf_slabs[1] += 1
                        continue
                    else:
                        weigth = types[lbl]/mass_slabs[1]
                weigths.append(weigth)
            weights = np.array(weigths, dtype=float)
            
        
        #Save position and times to calculate velocity at the next iteration
        elif step == 0:
            at_pos_1old = at_pos
            time_1old = time
            
        #Calculate velocity with velocity Verlet and write velocity files: 'vel.xyz', 'VCOM.dat'
        elif step > 0:
            time_vel = (time + time_2old)/2
            dt = (time - time_2old)
            at_vel = (at_pos - at_pos_2old)/dt
            at_vel = at_vel*phy.AonFS_to_MonS

            at_pos_2old = at_pos_1old
            time_2old = time_1old
            
            at_pos_1old = at_pos
            time_1old = time
            
            VCOM = np.dot(weights, at_vel[n_gf_slabs[0]:-n_gf_slabs[0]])
            
            #Write the .xyz file: 'Vel.xyz'
            header = ('ATOMIC_VELOCITIES (m/s) at time (fs) = %15.10f' %time_vel)
            data = (at_lbl, at_vel)
            ut.write_xyz(file = f_VEL, data = data, header = header)
            
            #Write the .dat file: 'VCOM.dat'
            f_VCOM.writelines('%15.10f %25.15f %25.15f %25.15f \n' %(time_vel, VCOM[0], VCOM[1], VCOM[2]))
        
        #Calculate RCOM           
        RCOM = np.dot(weights, at_pos[n_gf_slabs[0]:-n_gf_slabs[0]])
        
        #Write the .dat file: 'RCOM.dat'
        f_RCOM.writelines('%15.10f %25.15f %25.15f %25.15f \n' %(time, RCOM[0], RCOM[1], RCOM[2]))
        
        #Step counts total iteration number
        step += 1

f_VEL.close()
f_VCOM.close()
f_RCOM.close()