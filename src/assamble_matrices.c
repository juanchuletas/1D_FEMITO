#include <stdio.h>
#include <stdlib.h>
#include <data_structures.h>
void AssambleGlobalMatrices(int Ne,int order,double *link_mat,double *s_mat,double *k_mat,double *v_mat,double *eMatS,double *eMatK,double *eMatV,struct Element *e)
{
	int p,l,m;
	int nodes=Ne+1;
	p = order+1;
	double coeff;

	for(int ei=0; ei<Ne; ei++)
	{
		coeff = 0.5*e[ei].h;
		for(int i=0; i<p*p; i++)
		{
			eMatS[i] = coeff*eMatS[i];
			eMatK[i] = eMatK[i]/coeff;
			eMatV[i] = coeff*eMatV[i];
		}
			for(int nu=0; nu<p; nu++)
			{
				l=link_mat[p*ei+nu];
				for(int mu=0; mu<p; mu++)
				{
					m=link_mat[p*ei+mu];
					s_mat[l*nodes+m] +=  eMatS[p*nu+mu];
					k_mat[l*nodes+m] +=  eMatK[p*nu+mu];
					v_mat[l*nodes+m] +=  eMatV[p*nu+mu];
					//printf("%d %d\n",l*nodos+m,p*nu+mu);

				}
			}	
	}


	s_mat[0] = 1.0000;
	s_mat[1] = 0.0000;
	s_mat[nodes] = 0.0000;

	s_mat[(nodes-1)*nodes+(nodes-1)] = 1.0000;
 	s_mat[(nodes-2)*nodes+(nodes-1)] = 0.0000;
 	s_mat[(nodes-1)*nodes+(nodes-2)] = 0.0000;
	
	k_mat[0] = 1.0000;
	k_mat[1] = 0.0000;
	k_mat[nodes] = 0.0000;

	k_mat[(nodes-1)*nodes+(nodes-1)] = 1.0000;
 	k_mat[(nodes-2)*nodes+(nodes-1)] = 0.0000;
 	k_mat[(nodes-1)*nodes+(nodes-2)] = 0.0000;
	
	v_mat[0] = 1.0000;
	v_mat[1] = 0.0000;
	v_mat[nodes] = 0.0000;

	v_mat[(nodes-1)*nodes+(nodes-1)] = 1.0000;
 	v_mat[(nodes-2)*nodes+(nodes-1)] = 0.0000;
 	v_mat[(nodes-1)*nodes+(nodes-2)] = 0.0000;

}

