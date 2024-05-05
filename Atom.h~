#ifndef ATOM_H
#define ATOM_H
#include "stdheader.h"
#include "getboxdim.h"
#include "getljcoeff.h"


typedef struct {
        int    id;
	double  x;
	double  y;
	double  z;
	double  q;
        double  eps;
        double  sigma;
} Atom;



Atom **atominfo(FILE *fp1,Ljcoeff *lj, Box box,int natoms);
void free_atom(Atom **atom,int natoms);
#endif
