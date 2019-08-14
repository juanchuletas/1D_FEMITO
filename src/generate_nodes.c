#include <stdio.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
extern int AtomicNumber(char targer[3]);
double FroeseFischer(int i, double z,double r,double h)
{
	return (exp(-r + h*(double)i)/z);

}
void AtomicMesh(double *x,int Ne,double r0,double rN,char atom[3])
{
	double z = (double)AtomicNumber(atom);
	printf("Atomic Number: %lf\n",z);
	double ri,f_r1,f_r2;
	double rmf=5.0;
	double hmf=1.0/32;
	double rmin = FroeseFischer(0,z,rmf,hmf);
	double rmax = FroeseFischer(Ne,z,rmf,hmf);


	f_r1 = (rN-r0)/(rmax-rmin);
	f_r2 = (r0*rmax-rN*rmin)/(rmax-rmin);

	for(int i=0;i<Ne+1;i++)
	{
		ri = FroeseFischer(i,z,rmf,hmf);
		x[i] = f_r1*ri + f_r2;
		//printf("Vertex value x[%d] = %lf\n",i,x[i]);
	}


}
void GenerateNodes(double *x,int Ne, double r0, double rN, char mesh[180],char atom[3])
{
	double hx;
	if(strcmp(mesh,"Uniform")==0)
	{
		hx = (rN-(r0))/(double)(Ne);
		       
		for(int i=0; i<Ne+1; i++)
		{
			x[i] = r0 + hx*i;
			//printf("Vertex value x[%d] = %lf\n",i,x[i]);
        
		}
	}
	else if(strcmp(mesh,"Froese-Fischer")==0)
	{
		
		AtomicMesh(x,Ne,r0,rN,atom);

	}
}
