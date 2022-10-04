# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 13:13:34 2022

@author: simula
"""

# Numerical constants
k_load = 2.4275587028e7
mol = 6.02214179e23

# Energy-----------------------
ry_to_eV = 13.6056980659
J_to_eV = 1.602176565e-19
ha_to_eV = 27.211396
ha_to_ry = 2
kjmol_to_eVA = 0.0103642688
#-----------------------------

# Lenght
bohr_to_A = 0.529177208

#Velocity
AonFS_to_MonS = 1e5

# Force
rybohr_to_eVA = 25.71104309541616
habohr_to_eVA = 51.42208619083232

#Force on single atom via pressure
# F_at [Ry/Bohr] = (Area[A^2]*Pressure[GPa]/N_at)*0.000242756 
# Pressure[GPa] = (F_at [Ry/Bohr]*N_at)/(Area[A^2]*0.000242756)

# Pressure
eVA3_to_GPa = 160.21766208