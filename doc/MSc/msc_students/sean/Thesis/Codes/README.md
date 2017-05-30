# Project_Genesis
A CCDT code for infinite matter using closed shell system space.

Important packages used are:
- Eigen (version 3.3.3 currently included in project)
- MPI
- OpenMP (see branch MPI_no_OMP to run without)
- Other standard C++ packages. The code relies heavily on the STD library

Current up-to-date code lies in the "CCD"-folder and is compiled with:

path/to/CCD/ make

A typical run example is:

path/to/CCD/ mpirun -n 4 ./CCD "HEG" 14 5 1 1e-16 3

which will run the program for:

- "HEG" - Homogeneous electron gas, currently valid arguments are HEG and MP (for Minnesota potential)
- 14 particles
- 5 shells
- Wigner-Seitz radius of 1.0 femtometer
- Energy precision of 1e-16 (you'll probably be wanting to lower this for large systems using full CCDT)
- CCDT, the last integer can go form 0 to 3, where
	- 0 = CCD (all diagrams)
	- 1 = CCDT-1
	- 2 = CCDT-2
	- 3 = CCDT (all diagrams)

Lastly, yes the project name is silly.
