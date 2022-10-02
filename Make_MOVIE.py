# -*- coding: utf-8 -*-
"""
Created on Sat Sep 10 16:09:07 2022

@author: simula
"""
import numpy as np
import os

system = '50H_5GPa_300K'
skip_step = 10
apply_PBC = True
replica = [2,1,1]

path_IN = './OUTPUT/'+system+'/'
path_OUT = path_IN
if not os.path.isdir(path_OUT):
   os.mkdir(path_OUT)

file_pos = path_IN+'Pos.xyz'
file_info = path_IN+'info.dat'

file_out = path_OUT + 'Mov_every'+str(skip_step)+'_('+str(replica[0])+'x'\
    +str(replica[1])+'x'+str(replica[2])+').xyz'


with open(file_info, 'r') as f:
    for line in f:
        if 'Cell' in line:
            cell = ([next(f).split() for i in range(3)])
            
cell = np.array(cell, dtype=float)
PBC_vec = np.diag(cell)

f_out = open(file_out, 'w')
with open(file_pos, 'r') as f:
    step = -1
    for line in f:
        if 'ATOMIC_POSITIONS' in line:
            at_list = []
            step += 1
            
            if step%skip_step == 0:
                for line in f:
                    tokens = line.split()
                    if len(tokens) != 4: break
                
                    at_list.append(tokens)
                
                #Apply PBC
                pos = np.array([at[1:] for at in at_list], dtype=float)
                
                if apply_PBC:
                    pos = pos - (pos//PBC_vec)*PBC_vec
                
                #Apply replica
                pos_out = np.array([])
                at_repl = []
                for i, rep in enumerate(replica):
                    for j in range(rep):
                        pos_temp = pos + j*cell[i]
                        
                        if pos_out.size == 0:
                            pos_out = pos_temp
                        else:
                            pos_out = np.append(pos_out, pos_temp, axis=0)
                            n = pos_out[:,0].size
                        
                        for at in at_list:
                            at_repl.append(at[0])
                
                f_out.writelines('%d \n' %n)
                f_out.writelines('ATOMIC_POSITIONS \n')
                for i, at in enumerate(at_list):
                    f_out.writelines('%4s %28.12f %28.12f %28.12f \n' \
                                     %(at[0], pos_out[i,0], pos_out[i,1], pos_out[i,2]))
                
f_out.close()               












