# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 13:09:00 2022

@author: simula
"""
import numpy as np
import sys
import units
import UtilsLow as UtLow
import UtilsHigh as UtHig

#Command line managment section------------------------------------------------
#Verify number of expected parameters
if len(sys.argv) != 4:
    if len(sys.argv) > 4:
        raise ValueError('Too many arguments: 3 expected!')
        
    elif len(sys.argv) < 4:
        raise ValueError('Too few arguments: 3 expected!')

#Get simulation input and output files and output directory
IN_fname = str(sys.argv[1])
OUT_fname = str(sys.argv[2])
DATA_path = str(sys.argv[3])

#Verify existance of simulation files and create output directory if not present
UtLow.check_file_existance(IN_fname)
UtLow.check_file_existance(OUT_fname)
UtLow.create_dir(DATA_path)

#Extraction initialization section---------------------------------------------
#Read simulation input and initialize data for extraction
iprint, dt, celldm, nat_qm, ntyp, cell, types = UtHig.read_CP_input(IN_fname)

#Reference reading of simulation output to perform compatibility check between input and output
atoms_lbl_ref, atoms_pos_ref, line, n_line = UtHig.read_ref_output(OUT_fname)

#Total number of atoms
nat_ref = len(atoms_lbl_ref)

#Number of qm atoms
nat_qm_ref = nat_qm

#Number of gf atoms
nat_gf_ref = nat_ref - nat_qm_ref

#Compatibility check between simulation input and output
atoms_qm_lbl_ref = [lbl for lbl in atoms_lbl_ref if 'gf' not in lbl]
compatibility = UtHig.compatibility_check(atoms_qm_lbl_ref, types)
if not compatibility: 
    raise ValueError('Inconsistency between simulation input \
                     \' %s \'an output \' %s \'' %(IN_fname, OUT_fname))

#Get simulation output information for data extraction
mass_slabs, up_slab_index, idx_z_sort = UtHig.get_slab_info(atoms_lbl_ref, atoms_pos_ref, types) 

#Output files creation section-------------------------------------------------
#Create output files
INFO_fname = DATA_path+'info.txt'
POS_fname = DATA_path+'Pos.xyz'
FUP_fname = DATA_path+'Forces_gf_up.dat'
FDW_fname = DATA_path+'Forces_gf_dw.dat'
TEM_fname = DATA_path+'Temp.dat'

#Open output files to write them in the main loop
f_POS = open(POS_fname, 'w')
f_FOR_up = open(FUP_fname, 'w')
f_FOR_up.writelines('Total forces on upper GF atoms: TimeStep(fs) Fx Fy Fz(GPa) \n')
f_FOR_dw = open(FDW_fname, 'w')
f_FOR_dw.writelines('Total forces on lower GF atoms: TimeStep(fs) Fx Fy Fz(GPa) \n')
f_TEM = open(TEM_fname, 'w')    
f_TEM.writelines('Temperature of GF atoms: TimeStep(fs) T(Kelvin) \n')  

#Data extraction section-------------------------------------------------------
#Data extraction loop: read simulation output and write data to output files        
with open(OUT_fname, 'r') as f:
    
    step = -1
    
    #Keep track of read lines and line number while parsing the file
    line = 'Start'
    n_line = 0
    while line != '':
                
        #Reading file section--------------------------------------------------
                #Collect temperature and forces on GF atoms
                tags = ['kt_gf', 'upper[GPa]', 'lower[GPa]']
                tags_map, line, n_read = UtHig.read_tag(f, tags, limit=nat_ref)
                n_line += n_read
                #Corrupted field, limit reading, EOF->correct eof
                if 'ERROR' in line: 
                    
                    #Correct EOF
                    if 'EOF' in line: 
                        line = ''
                        break
                    break
                
                #'after' -> read, 'refresh' -> skip, 'EOF' -> stop, else -> incompatible file, stop
                line, n_read = UtLow.skip_lines(f, 1)
                n_line += n_read
                if 'after' in line: 
                            #Enumerate steps
                            step += 1
                            line, n_read = UtLow.skip_lines(f, 1)
                            n_line += n_read
                            if 'ATOMIC_POSITIONS' in line:
                                #Collect QM atom positions
                                qm_atoms, line, n_read = UtHig.read_list(f, num_line = nat_qm_ref, num_col = 4)
                                n_line += n_read
                                #EOF or Inconsistent reading
                                if 'ERROR' in line: break
                            
                                qm_atom_lbl = []
                                qm_atom_pos = []
                                for at in qm_atoms:
                                    qm_atom_lbl.append(at[0])
                                    qm_atom_pos.append(at[1:])
                                
                                qm_atom_lbl = np.array(qm_atom_lbl, dtype=str)
                                qm_atom_pos = np.array(qm_atom_pos, dtype=float)
                            
                            #EOF                               
                            elif 'ERROR' in line: break
                        
                            #ATOMIC_POSITION expected
                            else: 
                                line = 'ERROR = ATOMIC_POSITION expected'
                                break
                            
                            line, n_read = UtLow.skip_lines(f, nat_qm+4)
                            n_line += n_read
                            if 'GF_ATOM_POSITIONS' in line: 
                                #Collect GF atom positions
                                gf_atoms, line, n_read = UtHig.read_list(f, num_line = nat_gf_ref, num_col = 3)
                                n_line += n_read
                                if 'ERROR' in line: break
                            
                                gf_atom_lbl = []
                                gf_atom_pos = []
                                for at in gf_atoms:
                                    gf_atom_lbl.append('Cgf')
                                    gf_atom_pos.append(at)
                                
                                gf_atom_lbl = np.array(gf_atom_lbl, dtype=str)
                                gf_atom_pos = np.array(gf_atom_pos, dtype=float)
                            
                            #EOF
                            elif 'ERROR' in line: break
                            
                            #GF_ATOM_POSITIONS expected
                            else: 
                                line = 'ERROR = GF_ATOM_POSITIONS expected'
                                break
                
                elif 'refresh' in line: #Don't collect data if 'refresh'
                            line, n_read = UtLow.skip_lines(f, int(2*nat_ref + 7))
                            n_line += n_read
                            continue
                
                #EOF
                elif 'ERROR' in line: break
                
                #after or refresh expected
                else: 
                    line = 'ERROR = after or refresh expected'
                    break
            
        #Data checking section-----------------------------------------------
                #Get atom number to check consistency in reading
                nat_gf = gf_atom_lbl.size
                nat_tot = nat_qm + nat_gf                                
                #Create all atom positions and labels
                atom_pos = np.append(qm_atom_pos, gf_atom_pos, axis=0)
                atom_lbl = np.append(qm_atom_lbl, gf_atom_lbl, axis=0)
                
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
                atom_pos = atom_pos*units.bohr_to_A
                
        #Writing output section------------------------------------------------
                #Write extracted data to files
                if step%1000 == 0: print('%d steps processed' %step)
                
                #Write the .xyz file: Pos.xyz
                header = ('ATOMIC_POSITIONS (angstrom) at time (fs) = %15.10f' %time)
                data = (atom_lbl, atom_pos)
                UtHig.write_xyz(file = f_POS, data = data, header = header)
            
                #Write the .dat files:
                f_FOR_up.writelines('%15.10f %25.15f %25.15f %25.15f \n' %(time, F_gf_up[0], F_gf_up[1], F_gf_up[2]))
                f_FOR_dw.writelines('%15.10f %25.15f %25.15f %25.15f \n' %(time, F_gf_dw[0], F_gf_dw[1], F_gf_dw[2]))
                f_TEM.writelines('%15.10f %25.15f \n' %(time, T_gf))
                
        #Prepare for next iteration--------------------------------------------     
                #Skip irrelevant data
                line, n_read = UtLow.skip_lines(f, nat_gf_ref+1)
                n_line += n_read
                
                #EOF
                if 'ERROR' in line: break
            
        #----------------------------------------------------------------------
    
    #If one iteration has been performed succesfully, 'info.txt' can be written
    if step > 1: 
                #Convert dict to list in order to write it to file
                types = list(types.values())
                data_dt = dt*iprint
                data_step = step*iprint
                
                #Count the steps
                step += 1
                
                #Collect data to be written into info_data
                info_data = [cell, types, nat_tot, nat_qm, nat_gf, data_dt, dt, step, data_step, \
                                (up_slab_index, mass_slabs[0]), (nat_tot - up_slab_index, mass_slabs[1])]
                
                UtHig.write_info(file_path = INFO_fname, info_data = info_data)
  
    
    #Raise an error if the simulation output file is not as expected
    if 'ERROR' in line:
        f_POS.close()
        f_FOR_up.close()
        f_FOR_dw.close()
        f_TEM.close()
        
        if '_read_tag' in line:
            ErrorMSG = UtHig.errors_read_tag(OUT_fname, tags_map, line, n_line)
            raise ValueError(ErrorMSG)
            
        elif '_read_list' in line:
            ErrorMSG = UtHig.errors_read_list(OUT_fname, line, n_line)
            raise ValueError(ErrorMSG)
            
        elif 'EOF' in line:
            ErrorMSG = 'EOF in file \' %s \' at line %d' %(OUT_fname, n_line)
            raise ValueError(ErrorMSG)
            
        elif 'expected' in line:
            ErrorMSG = 'In file \' %s \' at line %d ' %(OUT_fname, line) + line.split('=')[1]
            raise ValueError(ErrorMSG)
        
        
    elif line == '':
        print('File \' %s \': successful read %d lines with %d iterations processed' %(OUT_fname, n_line, step))

f_POS.close()
f_FOR_up.close()
f_FOR_dw.close()
f_TEM.close()
