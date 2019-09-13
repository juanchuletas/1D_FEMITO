#include <stdio.h>
#include <stdlib.h>
#include "data_structures.h"
//extern double *ElementalOverlapMatrix(int order);
extern double nPolExtrapolation(double *x,double *y, int N,double target);
double *ComputeDensityR(double *wfn,struct Element *e,struct Vertex *n,int Ne, int order, int atomicN, int sign)
{
	//MODULE TO COMPUTE THE DENSITY MATRIX
	//
	int Occ=atomicN/2;
	int nodes = Ne*order+1;
	int fembasis = nodes-2;
	int r_nodes= nodes-2;
	//double *eMatS = ElementalOverlapMatrix(order);
	double *rho = (double *)malloc(sizeof(double)*nodes);
	double r_wfn[r_nodes*Occ];

	for(int orb=0; orb<Occ; orb++)
	{
		int k=0;
		for(int i=0; i<Ne; i++)
		{
			for(int j=0; j<order; j++)
			{
				r_wfn[orb*fembasis + k] = (wfn[orb*fembasis + k]*sign)/e[i].n[j+1].x;
				printf("RWFN[%d] = %lf\n",k,r_wfn[orb*fembasis + k]);
				k++;
			}
		}
	}
	int poly = 7;
	double x[poly],y[poly];
	double target;
	switch(order)
	{
		case 1:
			x[0]=e[1].n[0].x; y[0]=r_wfn[0];
			x[1]=e[2].n[0].x; y[1]=r_wfn[1];
       			x[2]=e[3].n[0].x; y[2]=r_wfn[2];
                        x[3]=e[4].n[0].x; y[3]=r_wfn[3];
                        x[4]=e[5].n[0].x; y[4]=r_wfn[4];
                        x[5]=e[6].n[0].x; y[5]=r_wfn[5];
                        x[6]=e[7].n[0].x; y[6]=r_wfn[6];
			target = nPolExtrapolation(x,y,poly,0.0000000);
			break;
		case 2:
			x[0]=e[0].n[1].x; y[0]=r_wfn[0];
                        x[1]=e[0].n[2].x; y[1]=r_wfn[1];
                        x[2]=e[1].n[1].x; y[2]=r_wfn[2];
                        x[3]=e[1].n[2].x; y[3]=r_wfn[3];
                        x[4]=e[2].n[1].x; y[4]=r_wfn[4];
                        x[5]=e[2].n[2].x; y[5]=r_wfn[5];
                        x[6]=e[3].n[1].x; y[6]=r_wfn[6];
			target = nPolExtrapolation(x,y,poly,0.0000000);
			break;
		case 3:
			x[0]=e[0].n[1].x; y[0]=r_wfn[0];
                        x[1]=e[0].n[2].x; y[1]=r_wfn[1];
                        x[2]=e[0].n[3].x; y[2]=r_wfn[2];
                        x[3]=e[1].n[1].x; y[3]=r_wfn[3];
                        x[4]=e[1].n[2].x; y[4]=r_wfn[4];
                        x[5]=e[1].n[3].x; y[5]=r_wfn[5];
                        x[6]=e[2].n[1].x; y[6]=r_wfn[6];
			target = nPolExtrapolation(x,y,poly,0.0000000);
			break;


	}
	rho[0] = 2.0*target*target;
	for(int i=0; i<fembasis; i++)
	{
		rho[i+1] = 0.0;
		for(int orb=0; orb<Occ; orb++)
		{
			rho[i+1] += r_wfn[i + orb*fembasis]*r_wfn[i + orb*fembasis];
			rho[i+1] = 2.0*rho[i+1];

		}
	}
	rho[nodes-1] = 0.f;


	/*printf("DENSITY MATRIX FOR %d OCCUPIED ORBITALS\n",Norb);
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
	printf("END OF DENS MATRIX FUNCTION\n");*/
	return rho;
	
}
