#include <stdio.h>
#include <stdlib.h>
#include "data_structures.h"
#include "Elemental_Matrices.h"
extern double *ElementalOverlapMatrix(int order);

void GetDens(double *rho,double *wfn,int nodes)
{
	int fembasis = nodes-2;
	int Occ = 1;
	for(int i=0; i<fembasis; i++)
        {
                rho[i] = 0.0;
                for(int orb=0; orb<Occ; orb++)
                {
                        rho[i] += wfn[i + orb*fembasis]*wfn[i + orb*fembasis];
                        rho[i] = 2.0*rho[i];
			printf("Density = %lf\n",rho[i]);

                }
        }

}


double IntegrateDens(double *wfn,struct Element *e,int *link_mat,int Ne, int order)
{
	int nodes = Ne*order+1;
	int fembasis = nodes-2;
	int p=order+1;
	double rho[fembasis];
	double f[fembasis];
	GetDens(rho,wfn,nodes);

	double coeff;
        double *eMatS = ElementalOverlapMatrix(order);
	double qtot = 0.0;
	for(int ei=0; ei<Ne; ei++)
        {
                coeff = 0.5*e[ei].h;
                for(int mu=0; mu<p; mu++)
                {
                        f[order*ei + mu] = 0.0;
                        for(int nu=0; nu<p; nu++)
                        {
                                int m = link_mat[p*ei+nu];
                                //f[order*ei + mu] += coeff*rho[m]*eMatS[mu*p+nu];
                                qtot += coeff*rho[m]*eMatS[mu*p+nu];
				//qtot = qtot + f[order*ei+mu];
                                //printf("f[%d] = %lf\n",ei+mu,f[ei+mu]);
                        }
                }
        }


	return qtot;

}
