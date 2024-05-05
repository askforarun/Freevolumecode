#ifndef CUBECENTRES_H
#define CUBECENTRES_H
#include "stdheader.h"
#include "getboxdim.h"
typedef struct{
double *xc;
double *yc;
double *zc;
double vol;
} Centres;

Centres *cubecentres(Box box, int Nx,int Ny,int Nz);
void free_cubecentres(Centres *centres);
#endif 
