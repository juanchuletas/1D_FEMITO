#include <stdio.h>
#include <stdlib.h>

void ComputeDensityMatrix(double *dens_mat,double *wfn,int Ne, int order, int atomicN)
{
	//MODULE TO COMPUTE THE DENSITY MATRIX
	//
	int Norb=atomicN/2;
	int nodes = Ne*order+1;
	int fembasis = nodes-2;
	printf("DENSITY MATRIX FOR %d OCCUPIED ORBITALS\n",Norb);
	for(int sigma=0; sigma<fembasis; sigma++)
	{
		for(int nu=0; nu<fembasis; nu++)
		{
			for(int orb=0; orb<Norb; orb++)
			{
				//rho = rho + wfn[sigma + orb*fembasis]*wfn[sigma + orb*fembasis];
				dens_mat[sigma*fembasis + nu] += wfn[sigma + orb*fembasis]*wfn[nu + orb*fembasis];
				//printf("P[%d][%d] = %lf\n",sigma,nu,dens_mat[sigma*fembasis + nu]);

			}
		}
	}
}
