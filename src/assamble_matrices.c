#include <stdio.h>
#include <stdlib.h>
#include "data_structures.h"
#include "Elemental_Matrices.h"
void FirstGradePol(int Ne,int order,int *link_mat,double *s_mat,double *k_mat,double *v_mat,double *eMatS,double *eMatK,double *eMatV,struct Element *e,struct Vertex *pot)
{
	int p,l,m;
	int nodes=(order*Ne) + 1;
	double coeff;
	p = order+1;
	double v[p];
			
	for(int ei=0; ei<Ne; ei++)
	{
		coeff = 0.5*e[ei].h;
		v[0] = e[ei].pot[0].x;
		v[1] = e[ei].pot[1].x;
		FirstGradeElementalMatrices(ei,order,coeff,link_mat,eMatS,eMatK,eMatV,v);
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


}
void SecondGradePol(int Ne,int order,int *link_mat,double *s_mat,double *k_mat,double *v_mat,double *eMatS,double *eMatK,double *eMatV,struct Element *e,struct Vertex *pot)
{
        int p,l,m;
	int nodes=(order*Ne) + 1;
        double coeff;
        p = order+1;
	double vij[p];

        for(int ei=0; ei<Ne; ei++)
        {
                coeff = 0.5*e[ei].h;
		vij[0] = e[ei].pot[0].x;
		vij[1] = e[ei].pot[1].x;
		vij[2] = e[ei].pot[2].x;
                SecondGradeElementalMatrices(ei,order,coeff,link_mat,eMatS,eMatK,eMatV,vij);
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


}
void AssambleGlobalMatrices(int Ne,int order,int *link_mat,double *s_mat,double *k_mat,double *v_mat,double *v,struct Element *e,struct Vertex *pot)
{
	int p,l,m;
	int nodes=(order*Ne) + 1;
	p = order+1;
	double *eMatS,*eMatK,*eMatV;
	eMatS = (double *)malloc(sizeof(double)*(p*p));
	eMatV = (double *)malloc(sizeof(double)*(p*p));
	eMatK = (double *)malloc(sizeof(double)*(p*p));

	switch(order)
	{
		case 1:
			FirstGradePol(Ne,order,link_mat,s_mat,k_mat,v_mat,eMatS,eMatK,eMatV,e,e->pot);
			break;
		case 2: 
			SecondGradePol(Ne,order,link_mat,s_mat,k_mat,v_mat,eMatS,eMatK,eMatV,e,e->pot);
			break;
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


	free(eMatS);
	free(eMatK);
	free(eMatV);

}

