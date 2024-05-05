#include "stdheader.h"
#include "Atom.h"
#include "sortatomscell.h"


Cell **sortatomscell(Atom **atom,int natoms,Box box,int Nx, int Ny, int Nz,double delx, double dely, double delz){

int i=0,j=0,k=0;
int l=0,m=0,n=0;
int count=0;
double x,y,z;
int cx=0;
int cy=0;
int cz=0;
int c=0,cp=0;
//int maxcells=Nx+(Ny*(Nx+2))+(Nz*(Nx+2)*(Ny+2));
int maxcells=(Nx+1)+((Ny+1)*(Nx+2))+((Nz+1)*(Nx+2)*(Ny+2));


Cell **head,*curr;

head=malloc(sizeof(Cell *)*(maxcells+1));

for(i=0;i<=maxcells;i++){
head[i]=NULL;
}


//printf("%lf %lf %lf\n",delx,dely,delz);

Cell *addfrontnode(Cell *head,Atom *atom)
{
Cell *curr=malloc(sizeof(Cell));
curr->atom=atom;
//printf("%lf\n",curr->atom->x);
curr->next=head;
return curr;
}

for(i=0;i<natoms;i++){
cx=ceil((atom[i]->x-box.xlo)/delx);
cy=ceil((atom[i]->y-box.ylo)/dely);
cz=ceil((atom[i]->z-box.zlo)/delz);
//printf("%d %d %d\n",Nx,Ny,Nz);

if (fabs(atom[i]->x-box.xlo) <= 1e-5 && cx==0)   {cx=1;}
if (fabs(atom[i]->x-box.xhi) <= 1e-5 && cx==Nx+1){cx=Nx;}
if (fabs(atom[i]->y-box.ylo) <= 1e-5 && cy==0)   {cy=1;}
if (fabs(atom[i]->y-box.yhi) <= 1e-5 && cy==Ny+1){cy=Ny;}
if (fabs(atom[i]->z-box.zlo) <= 1e-5 && cz==0)   {cz=1;}
if (fabs(atom[i]->z-box.zhi) <= 1e-5 && cz==Nz+1){cz=Nz;}

c=scalarcellindex(cx,cy,cz,Nx,Ny,Nz);

if(cx<1){printf("%d %d %d %d %d %lf %lf\n",c,cx,cy,cz,atom[i]->id,box.xlo,atom[i]->x);exit(1);}
else if(cx>Nx){printf("%d %d %d %d %d %lf %lf\n",c,cx,cy,cz,atom[i]->id,box.xhi,atom[i]->x);exit(1);}
if(cy<1){printf("%d %d %d %d %d %lf %lf\n",c,cx,cy,cz,atom[i]->id,box.ylo, atom[i]->y);}
else if(cy>Ny){printf("%d %d %d %d %d %lf %lf\n",c,cx,cy,cz,atom[i]->id,box.xhi,atom[i]->y);exit(1);}
if(cz<1){printf("%d %d %d %d %d %lf %lf\n",c,cx,cy,cz,atom[i]->id,box.zlo, atom[i]->z);exit(1);}
else if(cz>Nz){printf("%d %d %d %d %d %lf %lf\n",c,cx,cy,cz,atom[i]->id,box.zhi,atom[i]->z);exit(1);}
//printf("%d %d %d %d %d\n",atom[i]->id,c,cx,cy,cz);
//printf("%d\n",c);
curr=addfrontnode(head[c],atom[i]);
head[c]=curr;
}

/*for(i=1;i<=Nx;i++){
for(j=1;j<=Ny;j++){
for(k=1;k<=Nz;k++){
c=scalarcellindex(i,j,k,Nx,Ny,Nz);
while(head[c]){
count=count+1;
if(count==1)
{
x=head[c]->atom->x;
y=head[c]->atom->y;
z=head[c]->atom->z;
}
head[c]=head[c]->next;
}
//printf("%d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",c,count,i,j,k,x,y,z,box.xlo,box.xhi,box.ylo,box.yhi,box.zlo,box.zhi);
count=0;
}
}
}*/




return head;
}






