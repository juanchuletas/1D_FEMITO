#include <stdio.h>
#include <stdlib.h>
double *ComputeDensity(double *wfn,int Ne,int order,int atomicN)
{
	int nodes = Ne*order + 1;
	int fembasis = nodes-1;
	int Occ = atomicN/2;
	double *rho = (double *)malloc(sizeof(double)*fembasis);

	for(int i=0; i<fembasis; i++)
	{
		for(int orb=0; orb<Occ; orb++)
		{
			rho[i] = wfn[i + orb*fembasis]*wfn[i + orb*fembasis];
			rho[i] = 2.0*rho[i];
		}
	}
	return rho;
}
