# Freevolumecode
Calculation of free volume using hard sphere probe technique. This code was used in 10.1039/D1SM01825F and was written around 2018 while I was a pHD student. This uses a linked cell framework and is implemented using linked lists in C programming language. 
Unfortunately I did not have time to format the code but works. 


   
Innerworkings of the code:

The soratomcell.c is the crux of the code which sorts atoms in the cells. 
The main function is called widom.c but does not do any widom insertion. The entire free volume code was written as a part of free energy code that I wrote. 
The number of atoms in each cell can be queried by modifying sortatomcells.c
Each cell is identified by a scalar value (c in the code) and a vector, cx,cy,cz 
When looping, periodic boundary conditions/minimum image conventions are taken into account. 


Compile as 

gcc widom.c getnatoms.c getljcoeff.c  gettimestep.c getboxdim.c  Atom.c cubecentres.c particle.c delU.c sortatomscell.c min.c -o void

Use as 

./void  dump_300.lammpstrj allcoeffs_equalreactivity.txt 1.4 7

You can get the dumpfile from https://zenodo.org/records/11125570


The input for the code 
1) native lammps trajectories (i.e., dump files),
2) txt file containing atom diameters. The format of the text file is
   
   pair_coeff 1 1 lj/cut/coul/long 0.1094 3.39966950842 # c3

4) probe raidus and
5) the maximum LJ diameter of the atoms in the system (increasing this number will increase the number of neighbour cells scanned and hence will increase the compuational cost)

The output is 8 columns

1. number of neigbour cells scanned in each dimension (again can be increased when the last input is increased)
2. frame number (1 indexed)
3. FFV (fractional free volume) 
4. Free volume
5. Occupied volume
6. Total volume (the sum of 4+5)
7. Number of cubic free cells
8. Total number of cells

The code has the capability to output x,y,z coordinates of the centers of the free volume cubelets but I have supressed it. 

