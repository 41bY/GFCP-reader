# GFCP-reader

This project involves the creation of utilities to read and analyze the main output data file from coupled Green' s Function Car-Parrinello molecular dynamics (GFCP). The GFCP simulation software has been implemented as a plug-in within the Quantum Espresso suite, which is an open source suite to perform ab-initio atomistic simulation. GFCP is based on the standard CP software of Quantum Espresso (cp.x) and it is currently under further development in order to make it completely user interfaceable.

## GFCP in a nutshell

The GFCP package is meant for the atomistic simulation of material interfaces which are defined by two surfaces mating togheter. GFCP is a multi-scale simulation method because it employs the standard CP package in Quantum Espresso to perform ab-initio simulation of the interfacial atoms while it employs classical Green' s function approach to solve the dynamics of the bulk atoms. Ab-initio simulation relies on Density Functional Theory to solve the Schroedinger equation for the interfacial atoms. These kind of simulations are very accurate but also very computationally expansive so that only relatively small atomic systems can be practically simulated. 
Conversely classical molecular dynamic methods, as Green' s function, are less accurate but very fast from a computational point of view thus allowing to simulate large atomic systems. In this sense the linking of the two methods can give accurate results, because the interfacial chemistry is captured ab-initio while representing in a physical meaningful way the response of the bulk atoms, without effectively increasing the computational cost. 

## Usage

1) If the user wish to perform the extraction directly from GFCP simulation, then he/she can type from command line:

`python Extract.py path-to-GFCP-input/GFCP-input.in path-to-GFCP-output/GFCP-output.out path-to-OUTDIR/`

- path-to-GFCP-input/GFCP-input.in must be an existing file corresponding to the GFCP simulation input. This input has the format of a normal CP input file and in particular it must contain the tags: 'iprint', 'dt', 'celldm', 'nat', 'ntyp', 'cell' in order to be compatible with GFCP. For a complete description of CP input file [see](https://www.quantum-espresso.org/Doc/INPUT_CP.html). File name does not matter, the format does;

- path-to-GFCP-output/GFCP-output.out must be an existing file corresponding to the GFCP simulation output of the input 'GFCP-input'. File name does not matter, the format does;

- path-to-OUTDIR/ must be an output directory where the extracted data will be saved: 'Pos.xyz' 'Forces_gf_up.dat' 'Forces_gf_dw.dat' 'Temp.dat' will be stored. OUTDIR is created by the program if it is not found;

## Structure of the project

This project is divided as follow:
1) [Extract.py](https://github.com/41bY/GFCP-reader/blob/master/Extract.py)
2) [Read_input.py](https://github.com/41bY/GFCP-reader/blob/master/Read_input.py)
3) [Read_reference_output.py](https://github.com/41bY/GFCP-reader/blob/master/Read_reference_output.py)
4) [Read_output.py](https://github.com/41bY/GFCP-reader/blob/master/Read_output.py)
5) [Utils.py](https://github.com/41bY/GFCP-reader/blob/master/Utils.py)
6) [Testing](https://github.com/41bY/GFCP-reader/tree/master/Testing)
7) [Test_Read_input.py](https://github.com/41bY/GFCP-reader/blob/master/Test_Read_input.py)
8) [Test_Read_reference_output.py](https://github.com/41bY/GFCP-reader/blob/master/Test_Read_reference_output.py)
9) [Test_Read_output.py](https://github.com/41bY/GFCP-reader/blob/master/Test_Read_output.py)

1. [Extract.py](https://github.com/41bY/GFCP-reader/blob/master/Extract.py) has the purpose of extract the relevant simulation data directly from the GFCP simulation output file. These include the position of the atoms at each iteration 'Pos.xyz', the forces on the GF atoms at each iteration 'Force_gf_*.dat' and their temperature at each iteration 'Temp.dat'. The formats of the generated output files are standard file formats among the atomistic simulation community.
This program is not intended to perform heavy post-processing or data analysis on the GFCP output file but just to extract them and reorganize them in a more manageable form and standard form. This is done because the GFCP programs is somehow new and there are no available utilities to read and convert it in standard output data files. 
The way the program operates is to write the extracted data file while parsing the simulation output file. This is done because GFCP simulation output file is usually very large and it is not convenient or sometimes even possible to fully load it in memory. This program also perform 'safety checks' on the simulation output file while it is reading it to ensure the correctness of its structure. The main idea behind it being 'extract as much as possible' that is, stopping only when it cannot safely read the data, anyway generating a meaningful output (when possible).
This program relies on three 'main' modules [Read_input.py](https://github.com/41bY/GFCP-reader/blob/master/Read_input.py), [Read_reference_output.py](https://github.com/41bY/GFCP-reader/blob/master/Read_reference_output.py), [Read_reference_output.py](https://github.com/41bY/GFCP-reader/blob/master/Read_reference_output.py) togheter with a small utility module for trimming and unit conversion [Utils.py](https://github.com/41bY/GFCP-reader/blob/master/Utils.py). 

2. [Read_input.py](https://github.com/41bY/GFCP-reader/blob/master/Read_input.py) This module contains an head function 'read_CP_input' which controls the reading of the CP input file. This function relies on the other routines found inside this module and perform 'safety checks' to ensure the correctness of the CP input file format. [Test_Read_input.py](https://github.com/41bY/GFCP-reader/blob/master/Test_Read_input.py) is its respective testing module.

3. [Read_reference_output.py](https://github.com/41bY/GFCP-reader/blob/master/Read_reference_output.py) This module contains an head function 'read_ref_output' which controls the reading of the reference GFCP output file. Reference reading is necessary because GFCP output is usually a large file which is structured in iteration blocks. These iteration blocks have the same format that does not change within the GFCP output. The objective of this method is to read a reference block in order to get the indeces of each relevant tag within the block and thus speeding up the 'real' reading performed in [Extract.py](https://github.com/41bY/GFCP-reader/blob/master/Extract.py).
[Test_Read_reference_output.py](https://github.com/41bY/GFCP-reader/blob/master/Test_Read_reference_output.py) is its respective testing module. 

4. [Read_output.py](https://github.com/41bY/GFCP-reader/blob/master/Read_output.py) This module contains utility routines for the main program [Extract.py](https://github.com/41bY/GFCP-reader/blob/master/Extract.py), to efficiently read the GF-CP output file. These routines relies on parameters provided by the 'action' of the module [Read_reference_output.py](https://github.com/41bY/GFCP-reader/blob/master/Read_reference_output.py). [Test_Read_output.py](https://github.com/41bY/GFCP-reader/blob/master/Test_Read_output.py) is its respective testing module.

5. [Utils.py](https://github.com/41bY/GFCP-reader/blob/master/Utils.py) Small modules containing a trimming method and conversion factors for physical units. [Test_Utils.py](https://github.com/41bY/GFCP-reader/blob/master/Test_Utils.py) is its respective testing module.

6. [Testing](https://github.com/41bY/GFCP-reader/tree/master/Testing) This folder contains three folders for high level testing of the program and provide examples of output files.
In particular the folder [simulation_input](https://github.com/41bY/GFCP-reader/tree/master/Testing/simulation_input) contains a valid input file [cp.in](https://github.com/41bY/GFCP-reader/blob/master/Testing/simulation_input/cp.in) and various corrupted inputs 'cp-*'.
The folder [simulation_output](https://github.com/41bY/GFCP-reader/tree/master/Testing/simulation_output) contains a valid output file [gfcp.out](https://github.com/41bY/GFCP-reader/blob/master/Testing/simulation_input/cp.in) and various corrupted simulation outputs 'gfcp-*'.
The folder [OUTPUT](https://github.com/41bY/GFCP-reader/tree/master/Testing/OUTPUT) contains examples of the generated output files.

## Output data file formats

The output datafiles generated by [Extract.py](https://github.com/41bY/GFCP-reader/blob/master/Extract.py) are in agreement with standard atomistic simulation file formats, therefore making any post-processing easier with already existing post-processing programs.
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