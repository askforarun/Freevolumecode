# Freevolumecode
Calculation of free volume using hard sphere probe technique. This code was used in 10.1039/D1SM01825F 
This uses a linked cell framework and is implemented using linked lists in C.

The input for the code 
1) native lammps trajectories (i.e., dump files),
2) txt file containing atom diameters,
3) probe raidus and
4) the maximum LJ diamter of the atoms in the system (increasing this number will increase the number of cells scanned and hence will increase the compuational burden)
   
Unfortunately I did not have time to format the code. 
The soratomcell.c is the crux of the code which sorts atoms in the cells. 
The number of atoms in each cell can be queried by modifying sortatomcells.c
Each cell is designated a by a scalar and a vector (vector and mapped to a scalar and back) 

Compile as 

gcc widom.c getnatoms.c getljcoeff.c  gettimestep.c getboxdim.c  Atom.c cubecentres.c particle.c delU.c sortatomscell.c min.c -o void
Use as 

./void  dump_300.lammpstrj allcoeffs_equalreactivity.txt 0.7 4


The output is 8 columns

1. number of cells scanned in each dimension (again can be increased when the last input is increased)
2. frame number (1 indexed)
3. FFV (fractional free volume) 
4. Free volume
5. Occupied volume
6. Total volume
7. Number of cubic free cells
8. Total number of cells



