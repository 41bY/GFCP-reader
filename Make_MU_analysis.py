# -*- coding: utf-8 -*-
"""
Created on Sat Sep 10 16:09:07 2022

@author: simula
"""
import numpy as np
import os


system = '100H_1GPa_300K'
load = 5 #GPa

DATA_path = './OUTPUT/'+system+'/'
DATA_PES_path = DATA_path+'PES_analysis/'

file_gf_up = DATA_path+'Forces_gf_up.dat'
file_gf_dw = DATA_path+'Forces_gf_dw.dat'
file_Fpes = DATA_PES_path+'FonTraj_vs_time.dat'

KIN_file = DATA_path+'Kinetics.txt'
data_files = [file_gf_up, file_gf_dw, file_Fpes]

path_OUT = './OUTPUT/'+system+'/Force_spectra/'
if not os.path.isdir(path_OUT):
   os.mkdir(path_OUT)


#Calculate MU values and save them to file
f = open(KIN_file, 'a')
datas = []
f.write('\n')
for dataf in data_files:
    data = np.genfromtxt(dataf, skip_header=1)
    datas.append(data)
    
    F = np.average(data[:,1:3], axis=0)
    F_std = np.std(data[:,1:3], axis=0)
    
    #Distinguish between GF and PES forces (3 and 2 components respectively)
    if data[0].size < 4:
        Fz = load
        Fz_std = 0.0
        label = 'PES'
    else:
        Fz = np.average(data[:,3])
        Fz_std = np.std(data[:,3])
        label = dataf.split('/')[-1].split('Forces_')[1].split('.dat')[0]
        
    F = np.append(F, Fz)
    F_std = np.append(F_std, Fz_std)
    mu = np.abs(F[0]/F[2])
        
    f.write(label+'\t Fx Fy Fz(GPa) \t mu: \n')
    for i in range(3):
        f.write(' %10.5f +- %5.2f ' %(F[i], F_std[i]))
    f.write('\t %10.5f \n' %mu)
    
f.close()


#Calculate force spectra and save them to file
datas_fft = []
for data in datas:
    temp_data = np.fft.fft(data[:,1:], axis=0)
    dT = data[1,0] - data[0,0]
    temp_freq = np.fft.fftfreq(data[:,0].size, d=dT)
    ind = np.argwhere(temp_freq>0)
    x = temp_freq[ind][:,0]
    y = np.abs(temp_data[ind,:][:,0,:])**2
    x = np.array([x]).T
    fft_data = np.append(x, y, axis=1)
    datas_fft.append(fft_data)

for i, data in enumerate(datas_fft):
    np.savetxt(path_OUT+data_files[i].split('/')[-1],\
               data, header='Frequency(Phz) Forcespectralintensity(GPa^2)')