# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 13:09:00 2022

@author: simula
"""
import numpy as np
import os
from units import PhysConstants as phy

system = '100H_1GPa_300K'

# IN_path = str(sys.argv[1])
# OUT_path = str(sys.argv[2])

IN_path = '../../Data-CPGF-Diamond/'+system+'/'
OUT_path = '../../Data-CPGF-Diamond/'+system+'/GRNDIR/output/'
DATA_path = './OUTPUT/'+system+'/'

if not os.path.isdir(DATA_path.split(system)[0]):
   os.mkdir(DATA_path.split(system)[0])

if not os.path.isdir(DATA_path):
   os.mkdir(DATA_path)

fname_IN = IN_path + 'cp.in'
fname_GFOUT = OUT_path + 'gfcp.out'

fname_INFO = DATA_path+'info.txt'
fname_POS = DATA_path+'Pos.xyz'
fname_VEL = DATA_path+'Vel.xyz'
fname_FOR_up = DATA_path+'Forces_gf_up.dat'
fname_FOR_dw = DATA_path+'Forces_gf_dw.dat'
fname_TEM = DATA_path+'Temp.dat'

with open(fname_IN, 'r') as f:
    
    celldm = 1
    types = {}
    for line in f:
        
        if ('iprint' in line):
            token = line.split('=')[1]
            iprint = float(''.join([c for c in token if c in '0123456789']))           
        
        elif ('dt' in line): 
            token = line.split('=')[1]
            dt = float(''.join([c for c in token if c in '0123456789.']))
            dt = dt*2.4189*0.01 #fs
        
        elif ('celldm' in line): 
            token = line.split('=')[1]
            celldm = float(''.join([c for c in token if c in '0123456789.']))
            
        elif ('nat' in line): 
            token = line.split('=')[1]
            nat = int(''.join([c for c in token if c in '0123456789'])) 
            
        elif ('ntyp' in line): 
            token = line.split('=')[1]
            ntyp = int(''.join([c for c in token if c in '0123456789']))  
        
        elif ('CELL_PARAMETERS' in line):
            token = line.split()
            
            if len(token) > 1:
                unit = token[1]
                if unit == 'bohr':
                    cell_unit = phy.bohr_to_A
                elif unit == 'angstrom':
                    cell_unit = 1.0
                    
            else:
                cell_unit = celldm*phy.bohr_to_A
                
            cell = []
            for i in range(3):
                line = next(f)
                token = line.split()
                cell.append([float(c) for c in token])
            cell = cell_unit*np.array(cell)
            
        elif ('ATOMIC_SPECIES' in line):
            for i in range(ntyp):
                line = next(f)
                token = line.split()
                at = token[0]
                mass = float(token[1])
                types[at] = [mass, 0]

            
        elif ('ATOMIC_POSITIONS' in line):
            for i in range(nat):
                line = next(f)
                token = line.split()
                at = token[0]
                types[at][1] += 1

with open(fname_INFO, 'w') as f:
    
    f.writelines('Cell (angstrom):\n')
    for i in range(3):
        f.writelines('%15.10f %15.10f %15.10f' %(cell[i][0], cell[i][1], cell[i][2]))
        f.writelines('\n')

    f.writelines('\nAtom types (type|mass(a.m.u.)|#):\n')
    for ty in types.keys():
        f.writelines('%4s %10.4f %5d' %(ty, types[ty][0], types[ty][1]))
        f.writelines('\n')
        
    f.writelines('\nAtom number: %d \n' %nat)
    dt = iprint*dt
    f.writelines('\nData time step: %10.9f fs, simulation time step: %10.9f fs\n' %(dt, dt/iprint))
    

f_POS = open(fname_POS, 'w')
f_VEL = open(fname_VEL, 'w')

f_FOR_up = open(fname_FOR_up, 'w')
f_FOR_up.writelines('Total forces on upper GF atoms: TimeStep(fs) Fx Fy Fz(GPa) \n')

f_FOR_dw = open(fname_FOR_dw, 'w')
f_FOR_dw.writelines('Total forces on lower GF atoms: TimeStep(fs) Fx Fy Fz(GPa) \n')

f_TEM = open(fname_TEM, 'w')    
f_TEM.writelines('Temperature of GF atoms: TimeStep(fs) T(Kelvin) \n')      

pos_old = []    
with open(fname_GFOUT, 'r') as f:
    step = -1
    M_slabs = np.array([0.0, 0.0], dtype=float)
    order = True
  
    for line in f:
        write = False
        
        if ('kt_gf' in line): 
            token = line.split('=')[1].split()[0]
            T_gf = float(token)

        elif ('upper[GPa]' in line): 
            token = line.split('=')[1].split()
            F_gf_up = np.array([float(tk) for tk in token])

        elif ('lower[GPa]' in line): 
            token = line.split('=')[1].split()
            F_gf_dw = np.array([float(tk) for tk in token])
            
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
                
                
        # elif ('ATOMIC_VELOCITIES' in line):
        #     vel = []
        #     for i in range(nat):
        #         line = next(f)
        #         token = line.split()
        #         cords = [float(tk)*phy.HAunit_to_MperS for tk in token[1:]]
        #         vel.append([token[0]] + cords)
                
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


f_INFO = open(fname_INFO, 'a')
f_INFO.writelines('\nData step: %d, (simulation step is %dxData step) \n' %(step, int(iprint)))
f_INFO.writelines('\nSlabs info: (Slab \t Mass(a.m.u.))\n')
f_INFO.writelines('Lower \t %10.5f \n' %M_slabs[0])
f_INFO.writelines('Upper \t %10.5f \n' %M_slabs[1])

f_INFO.writelines('\nOrdering POS, VEL for QM atoms: QM upper slab start at %d \n'\
                  %i_split_QM)
for i in i_ord_QM:
    f_INFO.writelines('%4d ' %i)
    
f_INFO.writelines('\nOrdering POS, VEL for GF atoms: GF upper slab start at %d \n'\
                  %i_split_GF)
for i in i_ord_GF:
    f_INFO.writelines('%4d ' %i)

f_INFO.close()            
                
                
                
                
                
                
                
                
            
            
            
