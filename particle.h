#ifndef LJ_H
#define LJ_H
typedef struct {
double *x,*y,*z;
double natoms;
double *eps;
double *sigma;
}Particle;

Particle *createparticle(Particle *particle, double x, double y, double z);
Particle *particleinfo();
void free_particle(Particle *particle);
#endif
