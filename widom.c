#include "stdheader.h"
#include "getljcoeff.h"
#include "getboxdim.h"
#include "cubecentres.h"
#include "Atom.h"
#include "sortatomscell.h"
#include "particle.h"
#include <time.h>

//  gcc widom.c getnatoms.c getljcoeff.c  gettimestep.c getboxdim.c  Atom.c cubecentres.c particle.c delU.c sortatomscell.c min.c -lm -o void
 
//./a.out dump_0.45.lammpstrj paircoeff.txt 

//function declarations-----------

int getnatoms(FILE *fp);
double energy_diff(Box box, int natoms, Atom **atom, double cutoff,
		Particle *particle);
int gettimestep(FILE *fp);
double minimum(double x, double y,double z);
double diffenergy(Box box,Atom **atom, int natoms, double rc, double x, double y, double z,double rprobe);
//------------------------------------------------------------------------------

FILE *fp1, *fp2,*fp3;
int natoms = 0, timestep = 0, Nx =0, Ny =0, Nz =0, i = 0,j,k,counter = 0,
		counter_frame = 0;
int cx = 0, cy = 0, cz = 0, l, m, n, lp, mp, np,c,count,occ,occtemp;
double eps=0.0001;
double rc=3;
double delx,dely,delz,dx,dy,dz,dxs,dys,dzs;
double delUtemp = 0, exp_minusbetadelUtemp = 0, exp_minusbetadelUframe =0, exp_minusbetadelU = 0;
double delUtail,d,dtemp,rprobe=0,voidfrac;
double vol, void_vol;
double probevol,dprobe;
double maxsigma;
Particle *water;
int step=1;
int count=1;
Box box;
int N;
Atom **atom;
Centres *centre;
Atomparam *atomparam;
Ljcoeff *lj;
Cell **cell;
Cell *head;
char str1[100], str2[100];
char fil[100];
size_t nbytes = 100;
char *my_string;
int cbox;
ssize_t read;


int main(int argc, char *argv[]) {

int seed=time(NULL);
srand(seed);// Start the PRNG clock;
clock_t begin = clock();
        fp1 = fopen(argv[1], "r");
	fp2 = fopen(argv[2], "r");

        sscanf(argv[3], "%lf", &rprobe);
        sscanf(argv[4], "%lf", &maxsigma);
        
           
        dprobe = 2*rprobe;
        maxsigma =  maxsigma*0.5; 
	my_string = malloc(nbytes);
	natoms = getnatoms(fp1); //Get the number of atoms
//printf("%d\n",natoms);

	atomparam = getparameters(fp1, fp2, natoms);
  
      // printf("%lf\n",probevol);
	lj = Ljcoeffs(fp1, atomparam, natoms); //Get the Ljcoefficients epsilon and sigma for each atom in the system
	water = particleinfo(); // Get epsilon, sigma,Equilibrium position of the atoms,types of the atoms of the inserted particle-------------

	rewind(fp1); //Start from begining of the dump file
//read each line of the dump file
	while ((read = getline(&my_string, &nbytes, fp1))!= -1)  {

		
            
		sscanf(my_string, "%s %s", str1, str2);
                   
		if (sscanf(my_string, "%s %s", str1, str2) == 2) {
//------------------------------------------------------------------------------               
			if (strcmp(str1, "ITEM:") == 0 && strcmp(str2, "TIMESTEP") == 0) {
				counter = counter + 1;
				timestep = gettimestep(fp1);
                                
				//  printf("%d\n",timestep);
                             
			}
			if (strcmp(str1, "ITEM:") == 0 && strcmp(str2, "BOX") == 0) {
				box = boxdim(fp1);
				vol = (box.xhi - box.xlo) * (box.yhi - box.ylo)
						* (box.zhi - box.zlo);
                                 
				//printf("%lf\n",vol);
			}
//----------------Widom insertion begins here ------------------------------------------------------

			if (strcmp(str1, "ITEM:") == 0 && strcmp(str2, "ATOMS") == 0) {

if(count==counter){		
		
voidfrac=0;
occ=0;
//sprintf(fil, "fil_%d",counter);
//fp3=fopen(fil,"w+");
   
//printf("%lf\n",delUtail);
				atom = atominfo(fp1,lj, box,natoms);

//water=createparticle(water,1,1,4.3);//creates the particle at a desired location
//printf("%lf\n",water->lj->eps[0]);
//Atom sorting code should take atom and water info (ie coordinates) and sort them in boxes

                                     Nx = Nx(box.xhi,box.xlo,dprobe);
	                             Ny = Ny(box.yhi,box.ylo,dprobe);
	                             Nz = Nz(box.zhi,box.zlo,dprobe);

 	                             delx = delx(box.xhi,box.xlo,Nx);
	                             dely = dely(box.yhi,box.ylo,Ny);
	                             delz = delz(box.zhi, box.zlo,Nz);				


                                         N=ceil((maxsigma+rprobe)/dprobe);
                                       // printf("%d\n",N);
                           
                            
                           probevol=(double)(1.33)* 3.1416*rprobe*rprobe*rprobe;
//printf("%d %d %d %lf %lf %lf %lf %lf\n",Nx,Ny,Nz,delx,dely,delz,delx*dely*delz,probevol); 	

	centre = cubecentres(box, Nx, Ny, Nz);



//----------------------------------------------------------------------------------------------------------------------------
   //                            if(delx*dely*delz > probevol){
                                  
				//printf("%lf\n",delx*dely*delz);
				cell = sortatomscell(atom, natoms, box, Nx, Ny, Nz, delx, dely,delz);
			
				for (i = 0; i < Nx*Ny*Nz; i++) {
                                       //   printf("%d\n",i);
					cx = ceil((centre->xc[i] - box.xlo) / delx);
					cy = ceil((centre->yc[i] - box.ylo) / dely);
					cz = ceil((centre->zc[i] - box.zlo) / delz);
                                       
                                        cbox= scalarcellindex(cx, cy, cz, Nx, Ny, Nz);
                                         // printf("%d %d %d %d\n",cbox,cx,cy,cz);
                                    //  printf("%d %lf %lf %lf\n",cbox,centre->xc[i] ,centre->yc[i],centre->zc[i]);
                                    //  printf("%d %d %d %d %d %d\n",Nx,Ny,Nz,cx,cy,cz);
                                        // printf("%lf %lf %lf\n",c[i] ,centre->yc[i],centre->zc[i])
                                          
					for (l = cx - N; l <= cx + N; l++) {
						for (m = cy - N; m <= cy + N; m++) {
							for (n = cz - N; n <= cz + N; n++) {
                                                              //    printf("%d %d %d %d %d\n", i,c, lp, mp, np);
								if (l < 1) {lp = l + Nx;dxs = -(box.xhi-box.xlo);}
                                                                else if (l >Nx){lp=l-Nx;dxs = (box.xhi-box.xlo);} // Positive shift
                                                                else {lp=l;dxs=0;}
                                                                if (m < 1) {mp = m+ Ny;dys = -(box.yhi-box.ylo);}
                                                                else if (m>Ny){mp=m-Ny;dys = (box.yhi-box.ylo);}
                                                                else {mp=m;dys=0;}
                                                                if (n < 1) {np = n+ Nz;dzs = -(box.zhi-box.zlo);}
                                                                else if (n >Nz){np=n-Nz;dzs =(box.zhi-box.zlo);}
                                                                else {np=n;dzs=0;}
      								c= scalarcellindex(lp, mp, np, Nx, Ny, Nz);
                                                       //   printf("%d %d %d %d %d %d %d %d\n", i,c, l, m, n,lp,mp,np);          
       	                                                   if(cell[c]==NULL){dtemp=2;}
                                                               head=cell[c];
                                                               while(cell[c]){
                                                                dx=centre->xc[i]-(cell[c]->atom->x+dxs);
                                                                dy=centre->yc[i]-(cell[c]->atom->y+dys);
                                                                dz=centre->zc[i]-(cell[c]->atom->z+dzs);
                                                                 d=sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));
                                                              //  printf("%d %d %d %d %d %d %d %d\n", i,c, l, m, n,lp,mp,np);     
                                                               dtemp=d-(rprobe+0.5*(cell[c]->atom->sigma));

                                            
    //printf("%d %d %lf %lf %lf %lf %lf\n",c,cell[c]->atom->id,cell[c]->atom->x,cell[c]->atom->y,cell[c]->atom->z,d,dtemp); 
                                                 //   printf("%d %d %lf %lf %lf %lf\n",c,cell[c]->atom->id,cell[c]->atom->x,cell[c]->atom->y,cell[c]->atom->z,dtemp); 
                                                     // printf("%lf %lf\n",d,(rprobe+0.5*(cell[c]->atom->sigma)));
                                                              // printf("%lf %lf\n",d,(rprobe+0.5*(cell[c]->atom->sigma)));
                                                if(dtemp<eps){
                                                 // printf("%d %lf %lf %lf\n",i,cell[c]->atom->x,cell[c]->atom->y,cell[c]->atom->z);   
 
                                                              occ=occ+1; 
                                                  // printf("%d\n",i);
                                                        break;
                                                  }

                                                   cell[c]=cell[c]->next;
                                                              }                                                    
                                                             cell[c]=head;
                                                  if(dtemp<eps){
                                                     
                                                   //  printf("%lf\n",dtemp);
                                                        break;
                                                  }
                                                 	}                                               
                                                  if(dtemp<eps){
                                                      
                                             //        printf("%lf\n",d);
                                                        break;
                                                  }
						}
                                        if(dtemp<eps){
                                                     
                                               //      printf("%lf\n",d);
                                                       break;
                                                  }
					}
                                        
                                              if(dtemp>eps){
                                                    //rewind(fp3); 
                                                    //fprintf(fp3,"%lf %lf %lf\n",centre->xc[i] ,centre->yc[i] ,centre->zc[i]);
						    //printf("%lf %lf %lf %d %d\n",centre->xc[i] ,centre->yc[i] ,centre->zc[i],Nx*Ny*Nz-occ,counter);
                                              
                                                  }

                                  

			}
count = count+step;

voidfrac=(double) ((Nx*Ny*Nz)-occ)/(Nx*Ny*Nz);
printf("%d %d %lf %lf %lf %lf %d %d\n", N,counter,voidfrac,voidfrac*vol,vol-(voidfrac*vol),vol,Nx*Ny*Nz-occ,Nx*Ny*Nz);
}// timesteps
//------------------------------------------------------------------------------------------------*/

}
//}


//---------------------Widom insertion ends here---------------------

		}

	}
	// printf("%lf %d\n",-log(exp_minusbetadelU/counter)+2*ecor,counter);
//	free_getljcoeff(lj);
//	free(my_string);
//	free_getparameters(atomparam);
//	fclose(fp1);
//	fclose(fp2);

         clock_t end=clock();

double total_t = (double)(end - begin) / CLOCKS_PER_SEC;
printf("Total time taken by CPU: %lf\n", total_t);
	return 0;
}
