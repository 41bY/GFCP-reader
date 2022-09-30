# -*- coding: utf-8 -*-
"""
Created on Sat Sep 10 16:09:07 2022

@author: simula
"""
import numpy as np
from units import PhysConstants as phy

system = '100H_1GPa_300K'

DATA_path = './OUTPUT/'+system+'/'
INFO_file = DATA_path+'info.txt'
POS_file = DATA_path+'Pos.xyz'

COM_file = DATA_path+'COM.dat'
Z_file = DATA_path+'Slab_separation.dat'
VCOM_file = DATA_path+'VCOM.dat'
KIN_file = DATA_path+'Kinetics.txt'

n_at_plane = 16
types = ['C']
    
#Load atomic info
with open(INFO_file, 'r') as f:
    M = 0.0
    at = {}
    for line in f:
        if 'Atom types' in line:
            while(True):
                line = next(f)
                tokens = line.split()
                if len(tokens) != 3: break
                at[tokens[0]] = float(tokens[1])
                M += float(tokens[1])*int(tokens[2])
        
        elif 'Atom number' in line:
            nat = int(line.split(':')[1])

        elif 'Data time step' in line:
            token = line.split(':')[1].split('fs')[0]
            dt = float(''.join([c for c in token if c in '0123456789.']))
            
        elif 'Slabs info' in line:
            line = next(f)
            M_dw = float(line.split()[1])
            
            line = next(f)
            M_up = float(line.split()[1])             
            
        elif 'Ordering POS, VEL for QM atoms' in line:
            up_QM_init = int(line.split('slab start at')[1])
            tokens = next(f).split()
            order_QM = [int(i) for i in tokens]
            
        elif 'Ordering POS, VEL for GF atoms' in line:
            up_GF_init = int(line.split('slab start at')[1])
            tokens = next(f).split()
            order_GF = [int(i) for i in tokens]            
            
            
#Create relative COM, relative VEL, slab separation and save them to file
Rcom = []
Dz = []
Vcom = []
with open(POS_file, 'r') as f:
    step = -1
    for line in f:
        if 'ATOMIC_POSITIONS' in line:
            at_list = []
            for i in range(nat):
                line = next(f)
                tokens = line.split()
                if len(tokens) != 4: break
                if 'gf' in tokens[0]: break
                at_list.append(tokens)
                
            #Order atoms in ascending z-coordinate
            at_list_ord = np.array(at_list)[order_QM]
            wgts = [at[atom[0]] for atom in at_list_ord]
            
            #Calculate relative COM separation
            com_dw = np.average(at_list_ord[:up_QM_init,1:].astype(float), \
                                axis=0, weights=wgts[:up_QM_init])
            
            com_up = np.average(at_list_ord[up_QM_init:,1:].astype(float), \
                                axis=0, weights=wgts[up_QM_init:])
            com = com_up - com_dw
            Rcom.append(com)
            
            #Calculate slab separation
            z_pos_dw = np.array([float(atom[3]) for atom in at_list_ord[:up_QM_init]\
                                 if atom[0] in types])[-n_at_plane:]
                
            z_pos_up = np.array([float(atom[3]) for atom in at_list_ord[up_QM_init:]\
                                 if atom[0] in types])[:n_at_plane]                    
            
            z_dw = np.average(z_pos_dw)
            z_up = np.average(z_pos_up)
            Dz.append(z_up - z_dw)
            
            step += 1

Time_pos = dt*np.arange(0, step+1, 1)
Rcom = np.array(Rcom) - np.array([Rcom[0][0], Rcom[0][1], 0.0])
Zsep = np.array(Dz)
Vcom = ((Rcom[1:]-Rcom[:-1])/dt)*phy.AonFS_to_MonS #m/s 
Time_vel = (Time_pos[1:] + Time_pos[:-1])/2

avg_sep = np.average(Zsep)
std_sep = np.std(Zsep)

avg_Vcom = np.average(Vcom, axis=0)
std_Vcom = np.std(Vcom, axis=0)

with open(KIN_file, 'w') as f_out:
    f_out.writelines('Average slab (C-C) separation (Angstrom): %10.5f +- %5.2f \n' %(avg_sep, std_sep))
    f_out.writelines('Average COM relative velocity: Vx Vy Vz (m/s) \n')
    for i in range(3):
        f_out.writelines('\t %10.5f +- %5.2f' %(avg_Vcom[i], std_Vcom[i]))

np.savetxt(COM_file, np.column_stack((Time_pos, Rcom)),\
           header='Time(fs) Relative COM: X Y Z(Angstrom)')

np.savetxt(Z_file, np.column_stack((Time_pos, Zsep)),\
           header='Time(fs) C-C slab separation(Angstrom)')
    
np.savetxt(VCOM_file, np.column_stack((Time_vel, Vcom)),\
           header='Time(fs) Relative VCOM: Vx Vy Vz(m/s)')   

