# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 13:09:00 2022

@author: simula
"""
import numpy as np
import os
import sys
# from units import PhysConstants as phy
import Utils as ut


#Verify number of expected parameters
if len(sys.argv) != 4:
    
    if len(sys.argv) > 4:
        raise ValueError('Too many arguments')
        
    elif len(sys.argv) < 4:
        raise ValueError('Too few arguments')

IN_fname = str(sys.argv[1])
OUT_fname = str(sys.argv[2])
DATA_path = str(sys.argv[3])

ut.check_file_existance(IN_fname)
ut.check_file_existance(OUT_fname)
ut.create_dir(DATA_path)


INFO_fname = DATA_path+'info.txt'
POS_fname = DATA_path+'Pos.xyz'
VEL_fname = DATA_path+'Vel.xyz'
FUP_fname = DATA_path+'Forces_gf_up.dat'
FDW_fname = DATA_path+'Forces_gf_dw.dat'
TEM_fname = DATA_path+'Temp.dat'



input_data = ut.read_CP_input(IN_fname)
iprint = input_data[0]
dt = input_data[1]
celldm = input_data[2]
nat = input_data[3]
ntyp = input_data[4]
cell = input_data[5]
types = input_data[6]


with open(INFO_fname, 'w') as f:
    
    f.writelines('Cell (angstrom):\n')
    for i in range(3):
        f.writelines('%15.10f %15.10f %15.10f' %(cell[i][0], cell[i][1], cell[i][2]))
        f.writelines('\n')

    f.writelines('\nAtom types (type|mass(a.m.u.)|#):\n')
    for at in types:
        f.writelines('%4s %10.4f %5d' %(at[0], at[1], at[2]))
        f.writelines('\n')
        
    f.writelines('\nAtom number = %d \n' %nat)
    f.writelines('\nData time step = %10.9f fs\n' %(dt*iprint))
    f.writelines('\nSimulation time step = %10.9f fs\n' %(dt))
    

f_POS = open(POS_fname, 'w')
f_VEL = open(VEL_fname, 'w')

f_FOR_up = open(FUP_fname, 'w')
f_FOR_up.writelines('Total forces on upper GF atoms: TimeStep(fs) Fx Fy Fz(GPa) \n')

f_FOR_dw = open(FDW_fname, 'w')
f_FOR_dw.writelines('Total forces on lower GF atoms: TimeStep(fs) Fx Fy Fz(GPa) \n')

f_TEM = open(TEM_fname, 'w')    
f_TEM.writelines('Temperature of GF atoms: TimeStep(fs) T(Kelvin) \n')      

pos_old = []    
with open(OUT_fname, 'r') as f:
    step = -1
    M_slabs = np.array([0.0, 0.0], dtype=float)
    order = True
  
    for line in f:
        write = False
            
        tags = ['kt_gf', 'upper[GPa]', 'lower[GPa]']
        tags_map = ut.read_tag(f, tags)
        
        line = next(f)   
        if 'refresh' in line: #Don't collect data during refresh treatment
            for line in f:
                if 'GF_ATOM_VELOCITIES' in line:
                    break
            continue
        
        elif 'after' in line: #Collect data:
            line = next(f)
            
            if 'ATOMIC_POSITIONS' in line:
                atom_pos = ut.read_list(f, num_line = nat, check_consistency = True)
                atom_pos = np.array(atom_pos, dtype=str)
                
                if step < 0: #At first iteration, order the atoms following z-coordinate
                    z_pos = atom_pos[:,3].astype(dtype=float)
                    i_ord_QM = np.argsort(z_pos)
                    z_pos = z_pos[i_ord_QM]
                    dz = z_pos[1:] - z_pos[:-1]
                    i_split_QM = np.where(dz == dz.max())[0] + 1
                    
                    labels = np.array([row[0] for row in pos])
                    for i in range(nat):
                        if i < i_split_QM:
                            M_slabs[0] += types[labels[i_ord_QM[i]]][0]
                            
                        else:
                            M_slabs[1] += types[labels[i_ord_QM[i]]][0]
            
            
            
            else:
                raise ValueError('Line does not match the expectation: %s ,file is \
                                 corrupted!' %line)               
        
        else:
            raise ValueError('Line does not match the expectation: %s ,file is \
                             corrupted!' %line)
            
        
        T_gf = tags_map['kt_gf']
        T_gf = T_gf.split()[0]
        T_gf = float(T_gf)
        
        F_gf_up = tags_map['upper[GPa]']
        F_gf_up = F_gf_up.split()
        F_gf_up = np.array(F_gf_up, dtype=float) 
        
        F_gf_dw = tags_map['lower[GPa]']
        F_gf_dw = F_gf_dw.split()
        F_gf_dw = np.array(F_gf_dw, dtype=float)
            
        elif ('refresh' in line):
            for line in f:
                if 'GF_ATOM_VELOCITIES' in line:
                    break
            continue

        elif ('ATOMIC_POSITIONS' in line):
            pos = []
            for i in range(nat):
                line = next(f)
                token = line.split()
                cords = [float(tk)*phy.bohr_to_A for tk in token[1:]]
                pos.append([token[0]] + cords)
            
            if order:
                z_pos = np.array([row[3] for row in pos])
                i_ord_QM = np.argsort(z_pos)
                z_pos = z_pos[i_ord_QM]
                dz = z_pos[1:] - z_pos[:-1]
                i_split_QM = np.where(dz == dz.max())[0] + 1
                
                labels = np.array([row[0] for row in pos])
                for i in range(nat):
                    if i < i_split_QM:
                        M_slabs[0] += types[labels[i_ord_QM[i]]][0]
                        
                    else:
                        M_slabs[1] += types[labels[i_ord_QM[i]]][0]
                
                
        elif ('GF_ATOM_POSITIONS' in line):
            nat_gf = 0
            for i in range(nat):
                line = next(f)
                token = line.split()
                if(len(token) != 3): break
                nat_gf += 1
                cords = [float(tk)*phy.bohr_to_A for tk in token[:]]
                pos.append(['Cgf'] + cords)
                
            if order:
                z_pos = np.array([row[3] for row in pos[nat:]])
                i_ord_GF = np.argsort(z_pos)
                z_pos = z_pos[i_ord_GF]
                dz = z_pos[1:] - z_pos[:-1]
                i_split_GF = np.where(dz == dz.max())[0] + 1 + nat
                i_ord_GF = [it+nat for it in i_ord_GF]
                order = False
            
            #Calculate all atom 's velocities
            if len(pos_old) == len(pos):
                vel = []
                for i, at in enumerate(pos):
                    Vx = ((at[1] - pos_old[i][1])/dt)*phy.AonFS_to_MonS
                    Vy = ((at[2] - pos_old[i][2])/dt)*phy.AonFS_to_MonS
                    Vz = ((at[3] - pos_old[i][3])/dt)*phy.AonFS_to_MonS
                    vel.append([at[0], Vx, Vy, Vz]) #m/s
            
            pos_old = pos
            write = True
                
                
        # elif ('GF_ATOM_VELOCITIES' in line):
        #     for i in range(nat_gf):
        #         line = next(f)
        #         token = line.split()
        #         if(len(token) != 3): break
        #         cords = [float(tk)*phy.HAunit_to_MperS for tk in token[:]]
        #         vel.append(['Cgf'] + cords)
        #     write = True
            

            
        if write:
            
            step += 1
            if step%1000 == 0:
                print('%d steps processed' %step)
            
            
            f_POS.writelines('%d \n' %(nat+nat_gf))
            f_POS.writelines('ATOMIC_POSITIONS (angstrom) \n')
            for atom in pos:
                f_POS.writelines('%5s %23.15f %23.15f %23.15f \n' %(atom[0], atom[1], atom[2], atom[3]))

            
            if step > 0:
                f_VEL.writelines('%d \n' %(nat+nat_gf))
                f_VEL.writelines('ATOMIC_VELOCITIES (m/s) \n')
                for atom in vel:
                    f_VEL.writelines('%5s %23.15f %23.15f %23.15f \n' %(atom[0], atom[1], atom[2], atom[3]))            

            
            f_FOR_up.writelines('%14.11f %23.15f %23.15f %23.15f \n' %(step*dt, F_gf_up[0], F_gf_up[1], F_gf_up[2]))
            f_FOR_dw.writelines('%14.11f %23.15f %23.15f %23.15f \n' %(step*dt, F_gf_dw[0], F_gf_dw[1], F_gf_dw[2]))
            f_TEM.writelines('%14.11f %10.4f \n' %(step*dt, T_gf))
                
                       
f_POS.close()
f_VEL.close()
f_FOR_up.close()
f_FOR_dw.close()
f_TEM.close()       


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
                
                
                
                
                
                
                
                
            
            
            
