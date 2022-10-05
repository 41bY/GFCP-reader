# GFCP-reader

This project involves the creation of utilities to read and analyze the main output data file from coupled Green' s Function Car-Parrinello molecular dynamics (GFCP). The GFCP simulation software has been implemented as a plug-in within the Quantum Espresso suite, which is an open source suite to perform ab-initio atomistic simulation. GFCP is based on the standard CP software of Quantum Espresso (cp.x) and it is currently under further development in order to make it completely user interfaceable.

## GFCP in a nutshell

The GFCP package is meant for the atomistic simulation of material interfaces which are defined by two surfaces mating togheter. GFCP is a multi-scale simulation method because it employs the standard CP package in Quantum Espresso to perform ab-initio simulation of the interfacial atoms while it employs classical Green' s function approach to solve the dynamics of the bulk atoms. Ab-initio simulation relies on Density Functional Theory to solve the Schroedinger equation for the interfacial atoms. These kind of simulations are very accurate but also very computationally expansive so that only relatively small atomic systems can be practically simulated. 
Conversely classical molecular dynamic methods, as Green' s function, are less accurate but very fast from a computational point of view thus allowing to simulate large atomic systems. In this sense the linking of the two methods can give accurate results, because the interfacial chemistry is captured ab-initio while representing in a physical meaningful way the response of the bulk atoms, without effectively increasing the computational cost. 

## Usage

1) If the user wish to perform the extraction directly from GFCP simulation, then he/she can type from command line:

`python Extract_quantities.py path-to-GFCP-input/GFCP-input.in path-to-GFCP-output/GFCP-output.out path-to-OUTDIR/`

- path-to-GFCP-input/GFCP-input.in must be an existing file corresponding to the GFCP simulation input. This input has the format of a normal CP input file and in particular it must contain the tags: iprint, dt, celldm, nat, ntyp, cell, types in order to be compatible with GFCP. For a complete description of CP input file [see](https://www.quantum-espresso.org/Doc/INPUT_CP.html). File name does not matter, the format does;

- path-to-GFCP-output/GFCP-output.out must be an existing file corresponding to the GFCP simulation output of the input 'GFCP-input'. File name does not matter, the format does;

- path-to-OUTDIR/ must be an output directory where the extracted data will be saved: 'Pos.xyz' 'Forces_gf_up.dat' 'Forces_gf_dw.dat' 'Temp.dat' 'info.txt' will be stored. OUTDIR is created by the program if it is not found;

## Structure of the project

This project is divided as follow:
1) [Extract_quantities.py](https://github.com/41bY/GFCP-reader/blob/master/Extract_quantities.py)
2) [UtilsHigh.py](https://github.com/41bY/GFCP-reader/blob/master/UtilsHigh.py)
3) [UtilsLow.py](https://github.com/41bY/GFCP-reader/blob/master/UtilsLow.py)
4) [units.py](https://github.com/41bY/GFCP-reader/blob/master/units.py)
5) [Testing](https://github.com/41bY/GFCP-reader/tree/master/Testing)

1. [Extract_quantities.py](https://github.com/41bY/GFCP-reader/blob/master/Extract_quantities.py) has the purpose of extract the relevant simulation data directly from the GFCP simulation output file. These include the position of the atoms at each iteration 'Pos.xyz', the forces on the GF atoms at each iteration 'Force_gf_*.dat' and their temperature at each iteration 'Temp.dat'. Furthermore it will generate a small txt file named 'info.txt' which collects the main simulation parameters in a compact form. The formats of the generated output files are standard file formats among the atomistic simulation community, apart for the 'info.txt' which helps the user to have a quick check of the system.
This program is not intended to perform heavy post-processing or data analysis on the GFCP output file but just to extract them and reorganize them in a more manageable form and standard form. This is done because the GFCP programs is somehow new and there are no available utilities to read and convert it in standard output data files. 
The way the program operates is to write the extracted data file while parsing the simulation output file. This is done because GFCP simulation output file is usually very large and it is not convenient or sometimes even possible to fully load it in memory. This program also perform 'safety checks' on the simulation output file while it is reading it to ensure the correctness of its structure. The main idea behind it being 'extract as much as possible' that is, stopping only when it cannot safely read the data anyway generating a meaningful output (when possible).

2) [UtilsHigh.py](https://github.com/41bY/GFCP-reader/blob/master/UtilsHigh.py) This module contains useful 'high levels' functions which are used by the main program [Extract_quantities.py](https://github.com/41bY/GFCP-reader/blob/master/Extract_quantities.py). 
These functions are: reading and writing text based files in different ways and functions which have the purpose of making 'safety checks' and initialize some needed parameter in the main program [Extract_quantities.py](https://github.com/41bY/GFCP-reader/blob/master/Extract_quantities.py).

3) [UtilsLow.py](https://github.com/41bY/GFCP-reader/blob/master/UtilsLow.py) This module contains useful 'low levels' functions. With low levels we intend 'simple functions' which are prevalently called by more sophisticated functions inside [UtilsHigh.py](https://github.com/41bY/GFCP-reader/blob/master/UtilsHigh.py).

4) [units.py](https://github.com/41bY/GFCP-reader/blob/master/units.py) This module contains conversion factor among different physical units.

5) [Testing](https://github.com/41bY/GFCP-reader/tree/master/Testing) This folder contains various two type of testings. The first are executed on low level routines while the other ones are executed on high level routines. Because this second type of routines involve reading and writing text files, folder with ad-hoc modified simulation input and output are present in order to verify the correct behavior of the functions in the presence of corrupted or mismatched files.
In particular the folder [simulation_input](https://github.com/41bY/GFCP-reader/tree/master/Testing/simulation_input) contains a valid input file [cp.in](https://github.com/41bY/GFCP-reader/blob/master/Testing/simulation_input/cp.in) and various corrupted inputs 'cp-*'.
The folder [simulation_output](https://github.com/41bY/GFCP-reader/tree/master/Testing/simulation_output) contains a valid output file [gfcp.out](https://github.com/41bY/GFCP-reader/blob/master/Testing/simulation_input/cp.in) and various corrupted simulation outputs 'gfcp-*'.
The folder [OUTPUT](https://github.com/41bY/GFCP-reader/tree/master/Testing/OUTPUT) contains a example of the generated output files.

## Output data file formats

The output datafiles generated by [Extract_quantities.py](https://github.com/41bY/GFCP-reader/blob/master/Extract_quantities.py) are in agreement with standard atomistic simulation file formats, therefore making any post-processing easier with already existing post-processing programs.
The file format '.xyz' is structured as:

```
N
'Comment line, usually specifing units and iterations'
atom_1_label    x1  y1  z1
atom_2_label    x2  y2  z2
...             ..  ..  ..
atom_N_label    xN  yN  zN
.
.
[Repeated for every iteration]
.
.
```

The file format '.dat' is structured as:

```
'Comment line, usually specifing units and iterations'
f1           g1           h1           ...
f2           g2           h2           ...
..  ..  ..  ...
..  ..  ..  ...
f_last_step  g_last_step  h_last_step  ...

Where f, g, h, ... are functions whose value is defined at each simulation step.
```