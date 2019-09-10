#include <stdio.h>
#include <stdlib.h>
#include "data_structures.h"
extern double *ElementalHartreePotentialMatrix(int order,double *vh);
extern void FillZeroMat(double *matA,int N, int M);
void GetHartreeMatrix(int Ne,int order,int *link_mat,double *h_mat,struct Element *e,double *vh)
{
        int p,l,m;
        int nodes=(order*Ne) + 1;
        double coeff;
        p = order+1;
        double v[p];
        double *eMatHV;
        for(int ei=0; ei<Ne; ei++)
        {
                coeff = 0.5*e[ei].h;
		for(int i=0; i<p; i++)
		{
			v[i] = vh[link_mat[p*ei + i]];
		}
                eMatHV = ElementalHartreePotentialMatrix(order,v);
                for(int nu=0; nu<p; nu++)
                {
                        l=link_mat[p*ei+nu];
                        for(int mu=0; mu<p; mu++)
                        {
                                m=link_mat[p*ei+mu];
                                h_mat[l*nodes+m] +=  coeff*eMatHV[p*nu+mu];
                                //printf("%d %d\n",l*nodos+m,p*nu+mu);
                        }
                }
        }
        free(eMatHV);
}

void AssambleHartreePot(int Ne,int order,int *link_mat,double *vhij,struct Element *e,double *v)
{
	int nodes = Ne*order + 1;
	double *h_mat = (double *)malloc(sizeof(double)*nodes*nodes);
	FillZeroMat(h_mat,nodes,nodes);
	GetHartreeMatrix(Ne,order,link_mat,h_mat,e,v);
	/*switch(order)
	{
		case 1: FirstGradeHartree(Ne,order,link_mat,h_mat,e,v);
			break;
		case 2: //SecondGradeHartree(Ne,order,link_mat,h_mat,e,v);
			break;
		case 3: //SecondGradeHartree(Ne,order,link_mat,h_mat,e,v);
			break;
	}*/

	h_mat[0] = 1.0000;
        h_mat[1] = 0.0000;
        h_mat[nodes] = 0.0000;

        h_mat[(nodes-1)*nodes+(nodes-1)] = 1.0000;
        h_mat[(nodes-2)*nodes+(nodes-1)] = 0.0000;
        h_mat[(nodes-1)*nodes+(nodes-2)] = 0.0000;

	int l=0;
        int n = nodes;
        for(int i=1;i<n-1;i++)
        {
                for(int j=1;j<n-1;j++)
                {
                        vhij[l] = h_mat[i*n+j];
                        l++;
                }
	}


	free(h_mat);


}
