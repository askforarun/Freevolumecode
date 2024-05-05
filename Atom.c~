#include "stdheader.h"
#include "Atom.h"
#include "getboxdim.h"
Atom **atominfo(FILE *fp1,Ljcoeff *lj, Box box,int natoms) {
			
                                int i=0;
                                                         
                                 size_t nbytes = 100;
	                        char *my_string=NULL;
                         	Atom **atom=malloc(sizeof(Atom *)*natoms);
                           //    printf("%d\n",natoms);
                            for (i = 0; i < natoms; i++) {
                                        atom[i]= malloc(sizeof(Atom));
                                        
                                        atom[i]->eps=lj->eps[i];
                                        atom[i]->sigma=lj->sigma[i];      
					getline(&my_string, &nbytes, fp1);
					sscanf(my_string, "%d %*d %lf %lf %lf %lf\n",  &atom[i]->id, &atom[i]->x,
							&atom[i]->y, &atom[i]->z, &atom[i]->q);
//                     				printf("%d %d %lf %lf %lf %lf\n",i+1,atom[i]->id,atom[i]->x,atom[i]->y,atom[i]->z,atom[i]->sigma);
                                            if((atom[i]->x - box.xlo) < 0.00001){
                                         //printf("%lf %lf\n",atom[i]->x,box.xlo);
                                        atom[i]->x= atom[i]->x-box.xlo+box.xhi;
                                        // printf("%lf %lf\n",atom[i]->x,box.xhi);
                                          }
                                      if((atom[i]->x -  box.xhi) > 0.00001){
                                     //   printf("%lf %lf\n",atom[i]->x,box.xhi);
                                        atom[i]->x= atom[i]->x-box.xhi+box.xlo;
                                       //  printf("%lf %lf\n",atom[i]->x,box.xlo);
                                          }
                                

                                       if((atom[i]->y -  box.ylo) < 0.00001){
                                        // printf("%lf %lf\n",atom[i]->x,box.ylo);
                                            atom[i]->y= atom[i]->y-box.ylo+box.yhi;
                                        // printf("%lf %lf\n",atom[i]->x,box.yhi);
                                          }
                                        if((atom[i]->y -  box.yhi) > 0.00001){
                                        // printf("%lf %lf\n",atom[i]->x,box.xlo);
                                         atom[i]->y= atom[i]->y-box.yhi+box.ylo;
                                         //printf("%lf %lf\n",atom[i]->x,box.xlo);
                                          }
                                  
                                          if((atom[i]->z -  box.zlo) < 0.00001){
                                       //  printf("%lf %lf\n",atom[i]->x,box.xlo);
                                         atom[i]->z= atom[i]->z-box.zlo+box.zhi;
                                        // printf("%lf %lf\n",atom[i]->x,box.xlo);
                                          }
                                          if((atom[i]->z - box.zhi) > 0.00001){
                                        // printf("%lf %lf\n",atom[i]->x,box.xlo);
                                        atom[i]->z= atom[i]->z-box.zhi+box.zlo;
                                       //  printf("%lf %lf\n",atom[i]->x,box.xlo);
                                          }
                                    
				                              
//printf("%d %d %lf %lf %lf\n",i+1,atom[i]->id,atom[i]->x,atom[i]->y,atom[i]->z);
				}                              
				return atom;
                              
			
}




		
void free_atom(Atom **atom, int natoms)
{
int i;
for(i=0;i< natoms;i++){
free(atom[i]);
}
free(atom);
}
