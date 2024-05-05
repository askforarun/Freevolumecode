#include "stdheader.h"
#include "getljcoeff.h"
#include "getboxdim.h"
#include "cubecentres.h"
#include "Atom.h"
#include "sortatomscell.h"
#include "particle.h"
#include <time.h>

//  gcc readfile.c getnatoms.c getljcoeff.c  gettimestep.c getboxdim.c  Atom.c cubeatoms.c particle.c delU.c sortatomscell.c min.c -lm
 
//./a.out dump_0.45.lammpstrj paircoeff.txt 

//function declarations-----------

int getnatoms(FILE *fp, size_t nbytes, char *my_string);
double energy_diff(Box box, int natoms, Atom **atom, double cutoff,
		Particle *particle);
int gettimestep(FILE *fp, size_t nbytes, char *my_string);
double minimum(double x, double y,double z);
double diffenergy(Box box,Atom **atom, int natoms, double rc, double x, double y, double z,double rprobe);
//------------------------------------------------------------------------------

FILE *fp1, *fp2,*fp3,*fp4,*fp5;
int natoms = 0, timestep = 0, Nx =0, Ny =0, Nz =0, i = 0,j,k,counter = 0,
		counter_frame = 0,numvoids;
int cx = 0, cy = 0, cz = 0, l, m, n, lp, mp, np,c,count,occ,occtemp;
double eps=0.0001;
double rc=5;
double delx,dely,delz,dx,dy,dz,dxs,dys,dzs;
double delUtemp = 0, exp_minusbetadelUtemp = 0, exp_minusbetadelUframe =0, exp_minusbetadelU = 0;
double delUtail,d,dtemp,rprobe=0,voidfrac;
double vol, void_vol;
double probevol,dprobe;
double maxsigma;
Particle *water;
int step=1;
int count=1,countm=0;
Box box;
int N;
Atom **atom;
Atom **atom1;
Atomparam *atomparam;
Ljcoeff *lj;
Cell **cell;
Cell *head;
char str1[100], str2[100];
char fil[100];
size_t nbytes = 100;
char *my_string;
int cbox;
int *ind1, totindex1;
ssize_t read;

int main(int argc, char *argv[]) {

      fp1 = fopen(argv[1], "r");
	fp2 = fopen(argv[2], "r");

        sscanf(argv[3], "%lf", &rprobe);
        sscanf(argv[4], "%lf", &maxsigma);
        fp3 = fopen(argv[5], "r");

     


       

int seed=time(NULL);
srand(seed);// Start the PRNG clock;
clock_t begin = clock();
  
           
        dprobe = 2*rprobe;
        maxsigma =  maxsigma*0.5; 
	my_string = malloc(nbytes);
	
	natoms = getnatoms(fp1, nbytes, my_string); //Get the number of atoms

ind1=malloc(sizeof(int)*natoms);  
i=0; 

while ((read = getline(&my_string, &nbytes, fp3)) != -1) {
sscanf(my_string,"%d",&ind1[i]);
//printf("%d %d\n",i,ind1[i]);
i=i+1;
}
totindex1=i;

 
atomparam = getparameters(fp1, fp2, natoms);



	
 
 
	lj = Ljcoeffs(fp1, atomparam, natoms); //Get the Ljcoefficients epsilon and sigma for each atom in the system

	water = particleinfo(); // Get epsilon, sigma,Equilibrium position of the atoms,types of the atoms of the inserted particle-------------

	rewind(fp1); //Start from begining of the dump file

atom1=malloc(sizeof(Atom *)*natoms);



for (i = 0; i < natoms; i++) {
                                        atom1[i]= malloc(sizeof(Atom));
                                        }
                               
//read each line of the dump file

	while (!feof(fp1)) {

		getline(&my_string, &nbytes, fp1);
		sscanf(my_string, "%s %s", str1, str2);
		if (sscanf(my_string, "%s %s", str1, str2) == 2) {
//------------------------------------------------------------------------------               
			if (strcmp(str1, "ITEM:") == 0 && strcmp(str2, "TIMESTEP") == 0) {
				counter = counter + 1;
				timestep = gettimestep(fp1, nbytes, my_string);
				//  printf("%d\n",timestep);
			}
			if (strcmp(str1, "ITEM:") == 0 && strcmp(str2, "BOX") == 0) {
				box = boxdim(fp1, nbytes, my_string);
				vol = (box.xhi - box.xlo) * (box.yhi - box.ylo)
						* (box.zhi - box.zlo);
				//printf("%lf\n",vol);
                            
			}
//----------------Widom insertion begins here ------------------------------------------------------

			if (strcmp(str1, "ITEM:") == 0 && strcmp(str2, "ATOMS") == 0) {
     

if(count==counter){		
//if(counter==1){	
countm=0;
numvoids=0;
voidfrac=0;
occ=0;
sprintf(fil, "fil_%d",counter);
fp4=fopen(fil,"r");
i=0;
while ((read = getline(&my_string, &nbytes, fp4)) != -1) {

sscanf(my_string,"%lf %lf %lf",&atom1[i]->x,&atom1[i]->y,&atom1[i]->z);
//printf("%d %lf\n",counter,atom1[i]->x);
i=i+1;
}
numvoids=i;






				atom = atominfo(fp1, nbytes, my_string, lj, box,natoms);

//water=createparticle(water,1,1,4.3);//creates the particle at a desired location
//printf("%lf\n",water->lj->eps[0]);
//Atom sorting code should take atom and water info (ie coordinates) and sort them in boxes

                                     Nx = Nx(box.xhi,box.xlo,rc);
	                             Ny = Ny(box.yhi,box.ylo,rc);
	                             Nz = Nz(box.zhi,box.zlo,rc);

 	                             delx = delx(box.xhi,box.xlo,Nx);
	                             dely = dely(box.yhi,box.ylo,Ny);
	                             delz = delz(box.zhi, box.zlo,Nz);				


                                         N=ceil((maxsigma+rprobe)/dprobe);
                                       // printf("%d\n",N);
                           
                            
                           probevol=(double)(1.33)* 3.1416*rprobe*rprobe*rprobe;
//printf("%d %d %d %lf %lf %lf %lf %lf\n",Nx,Ny,Nz,delx,dely,delz,delx*dely*delz,probevol); 	

	//atom = cubeatoms(box, Nx, Ny, Nz);



//----------------------------------------------------------------------------------------------------------------------------
   //                            if(delx*dely*delz > probevol){
                                  
				//printf("%lf\n",delx*dely*delz);
				cell = sortatomscell(atom, natoms, box, Nx, Ny, Nz, delx, dely,delz);
			
				for (i = 0; i < numvoids; i++) {
                                       //   printf("%d\n",i);
					cx = ceil((atom1[i]->x - box.xlo) / delx);
					cy = ceil((atom1[i]->y - box.ylo) / dely);
					cz = ceil((atom1[i]->z - box.zlo) / delz);
                                       
                                        cbox= scalarcellindex(cx, cy, cz, Nx, Ny, Nz);
                                      //  printf("%d %d %d %d\n",cbox,cx,cy,cz);
                                    //  printf("%d %lf %lf %lf\n",cbox,atom->xc[i] ,atom->yc[i],atom->zc[i]);
                                    //  printf("%d %d %d %d %d %d\n",Nx,Ny,Nz,cx,cy,cz);
                                        // printf("%lf %lf %lf\n",c[i] ,atom->yc[i],atom->zc[i])
                                          
					for (l = cx - 1; l <= cx + 1; l++) {
						for (m = cy - 1; m <= cy + 1; m++) {
							for (n = cz - 1; n <= cz + 1; n++) {
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
       	                                               
                                                               head=cell[c];
                                                               while(cell[c]){
                                                                          
                                                                dx=atom1[i]->x-(cell[c]->atom->x+dxs);
                                                                dy=atom1[i]->y-(cell[c]->atom->y+dys);
                                                                dz=atom1[i]->z-(cell[c]->atom->z+dzs);
                                                                 d=sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));
                                                              //  printf("%d %d %d %d %d %d %d %d\n", i,c, l, m, n,lp,mp,np);     
                                                             for(k=0;k < totindex1; k++) {
                                                               if(d < rc  && cell[c]->atom->id==ind1[k]){
                                                                          //     printf("%d\n", ind1[k]);
                                                                   countm=countm+1;
                                                                    }
                                                                     }
    
                                                   cell[c]=cell[c]->next;
                                                              }                                                    
                                                             cell[c]=head;
                                                 
                                                 	}                                               
                                             
						}
                                     
					}
                                        
                                       
                                  

			}

sprintf(fil, "spec");
fp5=fopen(fil,"w+");
printf("%d %lf %lf\n",counter,(double) countm/(numvoids*1.33333*3.141593*rc*rc*rc), (double) totindex1/vol);
//fprintf(fp5,"%d %lf %lf\n",counter,(double) countm/(numvoids*1.33333*3.141593*rc*rc*rc), (double) totindex1/vol);
count = count+step;
voidfrac=(double) ((Nx*Ny*Nz)-occ)/(Nx*Ny*Nz);
//printf("%d %d %lf %lf %lf %lf %d %d\n", N,counter,voidfrac,voidfrac*vol,vol-(voidfrac*vol),vol,Nx*Ny*Nz-occ,Nx*Ny*Nz);
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
