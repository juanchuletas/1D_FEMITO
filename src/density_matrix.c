#include <stdio.h>
#include <stdlib.h>
//extern double *ElementalOverlapMatrix(int order);
void ComputeDensityMatrix(double *dens_mat,double *wfn,int Ne, int order, int atomicN)
{
	//MODULE TO COMPUTE THE DENSITY MATRIX
	//
	int Norb=atomicN/2;
	int nodes = Ne*order+1;
	int fembasis = nodes-2;
	int r_nodes= nodes-2;
	//double *eMatS = ElementalOverlapMatrix(order);
	double *rho = (double *)malloc(sizeof(double)*nodes);
	printf("DENSITY MATRIX FOR %d OCCUPIED ORBITALS\n",Norb);
	int k=0;
	for(int sigma=0; sigma<fembasis; sigma++)
	{
		for(int nu=0; nu<fembasis; nu++)
		{
			dens_mat[sigma*fembasis + nu] = 0.0;
			for(int orb=0; orb<Norb; orb++)
			{
				//rho = rho + wfn[sigma + orb*fembasis]*wfn[sigma + orb*fembasis];
				dens_mat[sigma*fembasis + nu] += wfn[sigma + orb*fembasis]*wfn[nu + orb*fembasis];
				//printf("P[%d] = %lf\n",sigma*fembasis + nu,dens_mat[sigma*fembasis + nu]);

			}
			//printf("P[%d] = %lf\n",sigma*fembasis + nu,dens_mat[sigma*fembasis + nu]);
		}
	}
	printf("END OF DENS MATRIX FUNCTION\n");
	
	rho[0] = 0.f;
	rho[nodes-1] = 0.f;
	printf("DENSITY AT NODES: \n");
	printf("rho[%d] = %lf\n",0,rho[0]);
	for(int i=1; i<=r_nodes; i++)
	{
		//rho[i]=0.0;
		for(int nu=0; nu<fembasis*fembasis; nu++)
		{
				rho[i] += dens_mat[(i-1) + nu*fembasis];
				//ho[i] = 2.0*rho[i];

		}
		printf("rho[%d] = %lf\n",i,rho[i]);
	}
	printf("rho[%d] = %lf\n",nodes-1,rho[nodes-1]);
	
}
