# GFCP-reader

This project involves the creation of utilities to read and analyze output data from coupled Green' s Function Car-Parrinello molecular dynamics (GFCP). The GFCP simulation software has been implemented as a plug-in within the Quantum Espresso suite, which is an open source suite to perform ab-initio atomistic simulation. GFCP is based on the standard CP software of Quantum Espresso and it is currently under further development in order to make it completely user interfaceable.

## GFCP in a nutshell

The GFCP package is meant for the atomistic simulation of interfaces which is defined by two surfaces mating togheter. GFCP is a multi-scale simulation method because it employs the standard CP package in Quantum Espresso to perform ab-initio simulation of the interfacial atoms while it employs classical Green' s function approach to solve the dynamics of the bulk atoms. Ab-initio simulation relies on Density Functional Theory to solve the Schroedinger equation for the interfacial atoms. These kind of simulations are very accurate but also very computationally expansive so that only relatively small atomic systems can be practically simulated. 
Conversely classical molecular dynamic methods, as Green' s function, are less accurate but very fast from a computational point of view thus allowing to simulate large atomic systems. In this sense the linking of the two methods can give accurate results, because the interfacial chemistry is captured ab-initio while representing in a physical meaningful way the response of the bulk atoms, without increasing to match the computational cost. 

## Structure of the project

This project is divided in three main programs:
1) [Extract_quantities.py](https://github.com/41bY/GFCP-reader/blob/master/Extract_quantities.py)
2) [Kinetic_analysis.py](https://github.com/41bY/GFCP-reader/blob/master/Kinetic_analysis.py)
3) [PES_analysis.py](https://github.com/41bY/GFCP-reader/blob/master/PES_analysis.py)

1. [Extract_quantities.py](https://github.com/41bY/GFCP-reader/blob/master/Extract_quantities.py) has the purpose of extract the relevant simulation data directly from the GFCP simulation output file gfcp.out. These include the position of the atoms at each iteration Pos.xyz, the forces on the GF atoms Force_gf.dat and the temperature of the system Temp.dat. Furthermore it will generate a small txt file named info.txt which collects the main simulation parameters in a compact form. This program is not intended to perform post-processing or analysis on the output data but just to extract and collect them in a more manageable form. The way the program operates is to write the extracted data file while parsing the simulation output file. This is done because GFCP simulation output file is usually very large and it is not convenient/possible to fully load it in memory. This program also perform 'safety checks' on the simulation output file while it is reading it to ensure the correctness of its structure. The main idea behind it being 'extract as much as possible' that is, stopping only when it cannot safely read the data anyway generating a meaningful output (when possible).

2. [Kinetic_analysis.py](https://github.com/41bY/GFCP-reader/blob/master/Kinetic_analysis.py) has the purpose of performing low level post processing on the data files extracted by [Extract_quantities.py](https://github.com/41bY/GFCP-reader/blob/master/Extract_quantities.py). In particular it reads 'Pos.xyz' and 'info.txt', perform post processing on them, and write 'RCOM.dat', 'Vel.xyz' and 'VCOM.dat'.

3. [PES_analysis.py](https://github.com/41bY/GFCP-reader/blob/master/PES_analysis.py) has the purpose of performing high level post processing on the data files 'info.dat', 'RCOM.dat' and on another data file 'grid.txt'. This extra data file is not provided by [Kinetic_analysis.py](https://github.com/41bY/GFCP-reader/blob/master/Kinetic_analysis.py) and [Extract_quantities.py](https://github.com/41bY/GFCP-reader/blob/master/Extract_quantities.py) because it is derived from another kind of atomistic simulation (constrained relaxation in particular) already implemented within the standard Quantum Espresso suite. It can be used independently from [Extract_quantities.py](https://github.com/41bY/GFCP-reader/blob/master/Extract_quantities.py) provided a correct 'Cell' is given in 'info.txt'.

The user can stop at any extraction/post-processing level depending on his/her objective. 
Furthermore the two post-processing programs [Kinetic_analysis.py](https://github.com/41bY/GFCP-reader/blob/master/Kinetic_analysis.py) and [PES_analysis.py](https://github.com/41bY/GFCP-reader/blob/master/PES_analysis.py) are meant for an independent usage from [Extract_quantities.py](https://github.com/41bY/GFCP-reader/blob/master/Extract_quantities.py). Their role is to provide post-processing on output datafiles which are compatible with the majority of available atomistic simulation softwares. This is made doable by the choice of the input data extensions. '.xyz' and '.dat' are well known data extensions within the atomistic simulation community, while 'grid.txt' and 'info.dat' can be easily generated.

## Usage

1) If the user wish to perform the extraction directly from GFCP simulation, then he/she can type from command line:
python Extract_quantities.py path-to-GFCP-input/GFCP-input.in path-to-GFCP-output/GFCP-output.out path-to-OUTDIR/

- path-to-GFCP-input/GFCP-input.in must be an existing file corresponding to the GFCP input. File name does not matter, the format does;

- path-to-GFCP-output/GFCP-output.out must be an existing file corresponding to the GFCP output derived from the simulation whose input is GFCP-input.in. File name does not matter, the format does;

- path-to-OUTDIR/ must be an output directory where 'Pos.xyz' 'Forces_gf_up.dat' 'Forces_gf_dw.dat' 'Temp.dat' 'info.txt' will be stored.OUTDIR is created by the program if it is not found;

2) If the user has already executed [Extract_quantities.py](https://github.com/41bY/GFCP-reader/blob/master/Extract_quantities.py) or he/she already has 'Pos.xyz' (maybe generated with other simulation softwares) 'info.txt' 'RCOM.dat' 'grid.txt', then he/she can use the post-processing programs. In particular:

- [Kinetic_analysis.py](https://github.com/41bY/GFCP-reader/blob/master/Kinetic_analysis.py) requires 'Pos.xyz' and 'info.txt' with 'Slabs info'.

- [PES_analysis.py](https://github.com/41bY/GFCP-reader/blob/master/PES_analysis.py) requires 'RCOM.dat' 'info.txt' with 'Cell' and 'grid.txt'.

Please see the structure of 'info.txt' to know the specific formats for the tags 'Slabs info' and 'Cell'.

To use [Kinetic_analysis.py](https://github.com/41bY/GFCP-reader/blob/master/Kinetic_analysis.py) type from command line:
python Kinetic_analysis.py path-to-info/info.txt path-to-Pos/Pos.xyz path-to-OUTDIR/

- path-to-info/info.txt must be an existing file containing the tag 'Slabs info' as specified in 'info.txt'. File name does not matter, the format does;

- path-to-Pos/Pos.xyz must be an existing file corresponding to the atomic coordinates of the simulation in agreement with the standard '.xyz'. File name does not matter, the format does;

- path-to-OUTDIR/ must be an output directory where 'Vel.xyz' 'RCOM.dat' 'VCOM.dat' will be stored. OUTDIR is created by the program if it is not found;

To use [PES_analysis.py](https://github.com/41bY/GFCP-reader/blob/master/PES_analysis.py) type from command line:
python PES_analysis.py path-to-info/info.txt path-to-grid/grid.txt path-to-RCOM/RCOM.dat path-to-OUTDIR/

- path-to-info/info.txt must be an existing file containing the tag 'Cell' as specified in 'info.txt'. File name does not matter, the format does;

- path-to-grid/grid.txt must be an existing file corresponding to the potential energy surface (PES) to be interpolated.
The format is: x y energy ,please see 'grid.txt'. File name does not matter, the format does;

- path-to-RCOM/RCOM.dat must be an existing file corresponding to the relative center of mass slab displacement which is generated by  [Kinetic_analysis.py](https://github.com/41bY/GFCP-reader/blob/master/Kinetic_analysis.py) or with other means.
The format is: time X Y Z ,please see 'RCOM.dat'. File name does not matter, the format does;

- path-to-OUTDIR/ must be an output directory where 'PES_DATA\' 'PES_Potential\' and 'PES_Forces\' will be stored. OUTDIR is created by the program if it is not found;