#ifndef LJCOEFF_H
#define LJCOEFF_H
#include "stdheader.h"
typedef struct {
double *eps;
double *sigma;
}Ljcoeff;

typedef struct {
double *eps;
double *sigma;
int maxtype;
int *type;
int *typecount;

}Atomparam;

Atomparam *getparameters(FILE *fp1,FILE *fp2,int natoms);
Ljcoeff *Ljcoeffs(FILE *, Atomparam *atomparam,int natoms);
void free_getljcoeff(Ljcoeff *lj);
void free_getparameters(Atomparam *atomparam);
#endif 



