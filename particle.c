#include "stdheader.h"
#include "particle.h"

Particle *particleinfo()
{
Particle *particle=(Particle *)malloc(sizeof(Particle));
particle->natoms=1;
particle->x=(double *)malloc(sizeof(double)*particle->natoms);
particle->y=(double *)malloc(sizeof(double)*particle->natoms);
particle->z=(double *)malloc(sizeof(double)*particle->natoms);
particle->eps=(double *)malloc(sizeof(double)*particle->natoms);
particle->sigma=(double *)malloc(sizeof(double)*particle->natoms);
particle->eps[0]=1;
particle->sigma[0]=1;
particle->x[0]=1;
particle->y[0]=1;
particle->z[0]=1;
return particle;
}



Particle *createparticle(Particle *particle, double x, double y, double z)
{
particle->x[0]=x;
particle->y[0]=y;
particle->z[0]=z;
//printf("%lf\n",particle->x[0]);
return particle;
}


void free_particle(Particle *particle)
{
free(particle->x);
free(particle->y);
free(particle->z);
free(particle->eps);
free(particle->sigma);
free(particle);
}
