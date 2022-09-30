# -*- coding: utf-8 -*-
"""
Created on Sat Sep 10 16:09:07 2022

@author: simula
"""
import numpy as np
from scipy.interpolate import Rbf
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import os

save = True
load = False
system = '100H_1GPa_300K'

GRID_file = '../Pes-processing/Physical-system/grid-'+system.split('H')[0]+'%.txt'
DATA_path = './OUTPUT/'+system+'/'
INFO_file = DATA_path+'info.txt'
COM_file = DATA_path+'COM.dat'

n_replicas = 4
density = 50
npoints = 100
h = 0.1
tol = 1e-5
maxstep = 99999


if load:
    outdir = DATA_path+'PES_analysis/'
    
    xSave = np.genfromtxt(outdir+'/grid_plot_x.dat',skip_header=1)
    ySave = np.genfromtxt(outdir+'/grid_plot_y.dat',skip_header=1)
    Energy_repl = np.genfromtxt(outdir+'/grid_plot_E.dat',skip_header=1)
    Rcom = np.genfromtxt(outdir+'/Traj.dat',skip_header=1)
    Vvst = np.genfromtxt(outdir+'/PESonTraj_vs_time.dat',skip_header=1)
    Fvst = np.genfromtxt(outdir+'/FonTraj_vs_time.dat',skip_header=1)
    MEP = np.genfromtxt(outdir+'/MEP.dat',skip_header=1)
    
    xGrid_repl, yGrid_repl = np.meshgrid(xSave,ySave)
    xGrid_repl = xGrid_repl.T
    yGrid_repl = yGrid_repl.T
    
    
    # x_begin = xSave[0]
    # x_end = xSave[-1]
    Time = Vvst[:,0]
    V = Vvst[:,1]
    Fx = Fvst[:,1]
    Fy = Fvst[:,2]  

else:
    
    #Compute PES---------------------------------------------------------------
    #Load grid and cell data
    grid = np.loadtxt(GRID_file)
    with open(INFO_file, 'r') as f:
        cell = []
        for line in f:
            if 'Cell' in line:
                for i in range(3):
                    line = next(f)
                    tokens = [float(token) for token in line.split()]
                    cell.append(tokens)
    
    cell = np.array(cell)
    cell[0,:] /= 4
    cell[1,:] /= 2
    
    #Fold grid using PBC of cell
    grid[:,:2] = grid[:,:2] - (grid[:,:2]//np.array([cell[0,0], cell[1,1]]))*np.array([cell[0,0], cell[1,1]])
    grid = np.unique(grid, axis=0)
    
    #Shift grid centering to maximum of PES
    emin = grid[:,2].min()
    emax = grid[:,2].max()
    XYmax = grid[np.argwhere(grid[:,2] == emax)[0,0],:2]
    grid[:,:2] -= XYmax + np.array([0.0, cell[1,1]/2])
    grid[:,2] -= emin
    
    #Replicates grid
    if n_replicas%2 == 0:
        nMin = int(1-n_replicas/2)
        nMax = int(1+n_replicas/2)
    else:
        nMin = int(-n_replicas/2)
        nMax = int(1+n_replicas/2)    
    
    rep_grid = np.array([])
    for i in range(nMin, nMax, 1):
        for j in range(nMin, nMax, 1):
            replica = grid + np.array([i*cell[0,0], j*cell[1,1], 0.0])
            if(rep_grid.size == 0): 
                rep_grid = replica
                continue
            rep_grid = np.append(rep_grid, replica, axis=0)
    
    #Construct interpolator (centered on Maximum)
    rbf = Rbf(rep_grid[:, 0], rep_grid[:, 1], rep_grid[:, 2], function='quintic')
    
    #Construct PES on denser grid (centered on Maximum)
    den = complex(density*1j)
    x_dense_grid, y_dense_grid = np.mgrid[0:cell[0,0]:den, 0:cell[1,1]:den]
    dense_grid = np.column_stack((x_dense_grid.flatten(), y_dense_grid.flatten()))
    Energy = rbf(dense_grid[:,0], dense_grid[:,1])
    PES = np.column_stack((dense_grid, Energy))
    
    #Construct MEP
    # iMin = np.argmin(Energy)
    # XYmin = dense_grid[iMin,:]
    # Start = np.array([1.2, 1.2])
    Start = np.array([0.0, 0.7])
    Stop = Start + np.array([cell[0,0], 0.0])
    
    g = np.linspace(0, 1, npoints)
    curve = np.array([(Stop[0]-Start[0])*g + Start[0], np.ones(npoints)*Start[1]]).T
    
    dr = (curve[1:] - curve[:-1])
    dr = np.insert(dr, 0, 0.0, axis=0)
    dmin = np.max(np.max(dr, axis=0))
    # dx = np.min(dr[:,0])
    # dy = np.min(dr[:,1])
    # dl = np.sum(dr, axis=1)**2
    # dl = np.sqrt(dr)
    ds = 0.01
    
    dx = np.array([ds, 0.0])
    dy = np.array([0.0, ds])    
    
    # lcurve = np.cumsum(dr)
    for n in range(maxstep):
        
        
        # PBC_curve = curve - (curve//np.array([cell[0,0], cell[1,1]]))*np.array([cell[0,0], cell[1,1]])
        # V = rbf(PBC_curve[:,0], PBC_curve[:,1])

        
        curve_x = curve + dx
        PBC_curve_x = curve_x - (curve_x//np.array([cell[0,0], cell[1,1]]))*np.array([cell[0,0], cell[1,1]])
        Vpx = rbf(PBC_curve_x[:,0], PBC_curve_x[:,1])
        
        curve_y = curve + dy
        PBC_curve_y = curve_y - (curve_y//np.array([cell[0,0], cell[1,1]]))*np.array([cell[0,0], cell[1,1]])
        Vpy = rbf(PBC_curve_y[:,0], PBC_curve_y[:,1])                
        

        curve_x = curve - dx
        PBC_curve_x = curve_x - (curve_x//np.array([cell[0,0], cell[1,1]]))*np.array([cell[0,0], cell[1,1]])
        Vmx = rbf(PBC_curve_x[:,0], PBC_curve_x[:,1])
        
        curve_y = curve - dy
        PBC_curve_y = curve_y - (curve_y//np.array([cell[0,0], cell[1,1]]))*np.array([cell[0,0], cell[1,1]])
        Vmy = rbf(PBC_curve_y[:,0], PBC_curve_y[:,1])         
        
        
        dVx = -(Vpx - Vmy)/(2*ds)
        dVy = -(Vpy - Vmy)/(2*ds)
        dV = np.column_stack((dVx, dVy))
        # dV[0,:] = np.array([0.0, 0.0])
        # dV[-1,:] = np.array([0.0, 0.0])
        
        curve_try = curve + h*dV
        # PBC_curve_try = curve_try - (curve_try//np.array([cell[0,0], cell[1,1]]))*np.array([cell[0,0], cell[1,1]])
        # V_try = rbf(PBC_curve_try[:,0], PBC_curve_try[:,1])
        
        curve_try_x = curve_try + dx
        PBC_curve_try_x = curve_try_x - (curve_try_x//np.array([cell[0,0], cell[1,1]]))*np.array([cell[0,0], cell[1,1]])
        Vpx_try = rbf(PBC_curve_try_x[:,0], PBC_curve_try_x[:,1])
        
        curve_try_y = curve_try + dy
        PBC_curve_try_y = curve_try_y - (curve_try_y//np.array([cell[0,0], cell[1,1]]))*np.array([cell[0,0], cell[1,1]])
        Vpy_try = rbf(PBC_curve_try_y[:,0], PBC_curve_try_y[:,1])
        
        curve_try_x = curve_try - dx
        PBC_curve_try_x = curve_try_x - (curve_try_x//np.array([cell[0,0], cell[1,1]]))*np.array([cell[0,0], cell[1,1]])
        Vmx_try = rbf(PBC_curve_try_x[:,0], PBC_curve_try_x[:,1])
        
        curve_try_y = curve_try - dy
        PBC_curve_try_y = curve_try_y - (curve_try_y//np.array([cell[0,0], cell[1,1]]))*np.array([cell[0,0], cell[1,1]])
        Vmy_try = rbf(PBC_curve_try_y[:,0], PBC_curve_try_y[:,1])        
        
        dVx_try = -(Vpx_try - Vmx_try)/(2*ds)
        dVy_try = -(Vpy_try - Vmy_try)/(2*ds)
        dV_try = np.column_stack((dVx_try, dVy_try))
        # dV_try[0,:] = np.array([0.0, 0.0])
        # dV_try[-1,:] = np.array([0.0, 0.0])
        
        curve_old = curve
        curve = curve + h*((dV+dV_try)/2)
        
        #Reparametrize curve to avoid clots
        dl = curve[1:] - curve[:-1]
        dl = np.insert(dl, 0, [0.0, 0.0], axis=0)
        dl = np.cumsum(np.sqrt(dl[:,0]**2 + dl[:,1]**2))
        dl = dl/dl[-1]
        
        x_int = interp1d(dl, curve[:,0])
        y_int = interp1d(dl, curve[:,1])
        curve = np.column_stack((x_int(g), y_int(g)))
        
        # curve[0] = Start
        # curve[-1] = Stop
        
        
        dist = curve - curve_old
        dist = dist[1:] - dist[:-1]
        dist = np.insert(dist, 0, [0.0, 0.0], axis=0)**2
        L = np.sum(np.sqrt(np.sum(dist, axis=1)))
        
        if L < tol:
            print('Tolerance threeshold reached in %d steps' %n)
            break
    
    n_rep = int(np.max(curve[:,0])//cell[0,0] + 1)
    s_PES = PES
    for i in range(1,n_rep+1):
        shifted_PES = PES + i*np.array([cell[0,0], 0.0, 0.0])
        s_PES = np.append(s_PES, shifted_PES, axis=0)
    
    #Trasform 1d -> 2d arrays to use contourf
    xGrid = np.reshape(s_PES[:,0], (density*(1+n_rep), density))
    yGrid = np.reshape(s_PES[:,1], (density*(1+n_rep), density))
    Energy = np.reshape(s_PES[:,2], (density*(1+n_rep), density))
    
    fig = plt.figure(figsize=(15, 5), dpi=300)
    plt.rcParams['xtick.labelsize'] = 18
    plt.rcParams['xtick.major.size'] = 5
    plt.rcParams['ytick.labelsize'] = 18
    plt.rcParams['ytick.major.size'] = 5
    ax = fig.add_subplot(111)
    ax.set_aspect('equal')
    
    print(Energy.max())
    plt.contourf(xGrid, yGrid, Energy, levels=np.linspace(-0.01, 0.2825, 300), cmap=plt.cm.RdYlBu_r)
    ax.plot(curve[:,0], curve[:,1], color='yellow', linestyle='dashed')
    
    DataCOM = np.loadtxt(COM_file, skiprows=1)
    Time = DataCOM[:,0]
    Rcom = DataCOM[:,1:3]      
    
    #Center Rcom as the PES
    Rcom += -Rcom[0,:] + np.array([0.0, cell[1,1]/2])
    
    Xmax = Rcom[:,0].max()
    Ymax = Rcom[:,1].max()
    
    x_repl = int(Xmax//cell[0,0] + 2)
    y_repl = int(Ymax//cell[1,1] + 2)
    
    # x_repl = 5
    # Rcom = Rcom[Rcom[:,0] < x_repl*cell[0,0]]
    
    s_PES = PES
    MEP = curve
    
    for n in range(1,x_repl):
        PES_temp = PES + n*np.array([cell[0,0], 0.0, 0.0])
        s_PES = np.vstack((s_PES, PES_temp))
        
        MEP_temp = curve + n*np.array([cell[0,0], 0.0])
        MEP = np.vstack((MEP, MEP_temp))
    
    #Reparametrize MEP to smoothen it
    lmep = MEP[1:] - MEP[:-1]
    lmep = np.insert(lmep, 0, [0.0, 0.0], axis=0)
    lmep = np.cumsum(np.sqrt(lmep[:,0]**2 + lmep[:,1]**2))
    lmep = lmep/lmep[-1]
    
    x_int = interp1d(lmep, MEP[:,0])
    y_int = interp1d(lmep, MEP[:,1])
    MEP = np.column_stack((x_int(g), y_int(g)))
    
    xGrid_repl = np.reshape(s_PES[:,0], (density*x_repl, density))
    yGrid_repl = np.reshape(s_PES[:,1], (density*x_repl, density))
    Energy_repl = np.reshape(s_PES[:,2], (density*x_repl, density))
        
    #Construct potential and force on trajectory
    PBC_Rcom = Rcom - (Rcom//np.array([cell[0,0], cell[1,1]]))*np.array([cell[0,0], cell[1,1]])
    V = rbf(PBC_Rcom[:,0], PBC_Rcom[:,1])
    
    dr = 0.001*max([np.max(np.abs(Rcom[1:,0] - Rcom[:-1,0])), np.max(np.abs(Rcom[1:,1] - Rcom[:-1,1]))])
    Rcom_x = Rcom + np.array([dr, 0.0])
    Rcom_y = Rcom + np.array([0.0, dr])
    PBC_Rcom_x = Rcom_x - (Rcom_x//np.array([cell[0,0], cell[1,1]]))*np.array([cell[0,0], cell[1,1]])
    PBC_Rcom_y = Rcom_y - (Rcom_y//np.array([cell[0,0], cell[1,1]]))*np.array([cell[0,0], cell[1,1]])
    dVx = rbf(PBC_Rcom_x[:,0], PBC_Rcom_x[:,1])
    dVy = rbf(PBC_Rcom_y[:,0], PBC_Rcom_y[:,1])
    
    Fx = -((dVx-V)/dr)*10 #GPa
    Fy = -((dVy-V)/dr)*10
        
    #Save output files
    if save:
        outdir = DATA_path+'PES_analysis/'
        if not os.path.isdir(outdir):
            os.mkdir(outdir)
           
        xSave = np.unique(xGrid_repl, axis=1).flatten()
        ySave = np.unique(yGrid_repl, axis=0).flatten()
        np.savetxt(outdir+'/grid_plot_x.dat', xSave, header='Grid x component(Angstrom)')
        np.savetxt(outdir+'/grid_plot_y.dat', ySave, header='Grid y component(Angstrom)')
        np.savetxt(outdir+'/grid_plot_E.dat', Energy_repl, header='Grid energy: E(x,y) (J/m^2)')
        np.savetxt(outdir+'/Traj.dat', np.column_stack((Rcom[:,0], Rcom[:,1])), header='Centered CoM xy trajectory(Angstrom)')
        np.savetxt(outdir+'/PESonTraj_vs_time.dat', np.column_stack((Time, V)), header='Time(fs) \t V(J/m^2)')
        np.savetxt(outdir+'/FonTraj_vs_time.dat', np.column_stack((Time, Fx, Fy)), header='Time(fs) \t Fx \t Fy(GPa)')
        np.savetxt(outdir+'/MEP.dat', MEP, header='MEP X Y(Angstrom)')
        
    
#Plot Rcom on the PES    
fig = plt.figure(figsize=(15, 5), dpi=300)
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['xtick.major.size'] = 5
plt.rcParams['ytick.labelsize'] = 18
plt.rcParams['ytick.major.size'] = 5
ax = fig.add_subplot(111)
ax.set_aspect('equal')
    
plt.contourf(xGrid_repl, yGrid_repl, Energy_repl, levels=np.linspace(-0.01, 0.2825, 300), cmap=plt.cm.RdYlBu_r)
ax.plot(Rcom[:,0], Rcom[:,1], color='black', ls='solid', lw=0.5)
ax.plot(MEP[:,0], MEP[:,1], color='yellow', ls='dashed', lw=0.7)
plt.show()

plt.plot(Time, V)
plt.show()  

plt.plot(Time, Fx)
plt.show()    
            
plt.plot(Time, Fy)
plt.show()   
    
    


