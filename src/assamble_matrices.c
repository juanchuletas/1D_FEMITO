#include <stdio.h>
#include <stdlib.h>
#include "data_structures.h"
#include "Elemental_Matrices.h"

extern void ScalarXMatrix(double coeff,double *matA,double *matB,int N,int M);
extern void print_matrix(char *name,int m,int n,double *matA);
extern double *ElementalOverlapMatrix(int order);
extern double *ElementalKinectMatrix(int order);
extern double *ElementalPotentialMatrix(int order,double *v);
extern void ReduceMatrices(double *sij,double *kij,double *vij,double *uij,double *lij,double *s_mat,double *k_mat,double *v_mat,double *u_mat,double *l_mat,int nodes);
extern void FillZeroMat(double *A,int M,int N);
void ApplyBoundaryCond(double *s_mat,double *k_mat,double *v_mat,double *u_mat,double *l_mat,int nodes)
{

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


	u_mat[(nodes-1)*nodes+(nodes-1)] = 1.0000;
 	u_mat[(nodes-2)*nodes+(nodes-1)] = 0.0000;
 	u_mat[(nodes-1)*nodes+(nodes-2)] = 0.0000;

	l_mat[(nodes-1)*nodes+(nodes-1)] = 1.0000;
 	l_mat[(nodes-2)*nodes+(nodes-1)] = 0.0000;
 	l_mat[(nodes-1)*nodes+(nodes-2)] = 0.0000;
}
void FirstGradePol(int Ne,int order,int *link_mat,double *s_mat,double *k_mat,double *v_mat,struct Element *e,struct Vertex *pot)
{
	int p,l,m;
	int nodes=(order*Ne) + 1;
	double coeff;
	p = order+1;
	double v[p];
	double *eMatS,*eMatK,*eMatV;
	eMatS = ElementalOverlapMatrix(order);
        eMatK = ElementalKinectMatrix(order);
			
	for(int ei=0; ei<Ne; ei++)
	{
		coeff = 0.5*e[ei].h;
		v[0] = e[ei].pot[0].x;
		v[1] = e[ei].pot[1].x;
		eMatV = ElementalPotentialMatrix(order,v);
		for(int nu=0; nu<p; nu++)
		{
			l=link_mat[p*ei+nu];
			for(int mu=0; mu<p; mu++)
			{
				m=link_mat[p*ei+mu];
				s_mat[l*nodes+m] +=  coeff*eMatS[p*nu+mu];
				k_mat[l*nodes+m] +=  eMatK[p*nu+mu]/coeff;
				v_mat[l*nodes+m] +=  coeff*eMatV[p*nu+mu];
				//printf("%d %d\n",l*nodos+m,p*nu+mu);
			}
		}
	}
	printf("POT INTEGRAL = %lf\n",v_mat[0]);
	free(eMatV);
	free(eMatK);
	free(eMatS);


}
void SecondGradePol(int Ne,int order,int *link_mat,double *s_mat,double *k_mat,double *v_mat,struct Element *e,struct Vertex *pot)
{
        int p,l,m;
	int nodes=(order*Ne) + 1;
        double coeff;
        p = order+1;
	double vij[p];
	double *eMatS,*eMatK,*eMatV;
	eMatS = ElementalOverlapMatrix(order);
        eMatK = ElementalKinectMatrix(order);

        for(int ei=0; ei<Ne; ei++)
        {
                coeff = 0.5*e[ei].h;
		vij[0] = e[ei].pot[0].x;
		vij[1] = e[ei].pot[1].x;
		vij[2] = e[ei].pot[2].x;
		eMatV = ElementalPotentialMatrix(order,vij);
                for(int nu=0; nu<p; nu++)
                {
                        l=link_mat[p*ei+nu];
                        for(int mu=0; mu<p; mu++)
                        {
                                m=link_mat[p*ei+mu];
				s_mat[l*nodes+m] +=  coeff*eMatS[p*nu+mu];
				k_mat[l*nodes+m] +=  eMatK[p*nu+mu]/coeff;
				v_mat[l*nodes+m] +=  coeff*eMatV[p*nu+mu];
                                //printf("%d %d\n",l*nodos+m,p*nu+mu);
                        }
                }
        }
	free(eMatV);
	free(eMatK);
	free(eMatS);


}
void ThirdGradePol(int Ne,int order,int *link_mat,double *s_mat,double *k_mat,double *v_mat,struct Element *e,struct Vertex *pot)
{
        int p,l,m;
	int nodes=(order*Ne) + 1;
        double coeff;
        p = order+1;
	double vij[p];
	double *eMatS,*eMatK,*eMatV;
	eMatS = ElementalOverlapMatrix(order);
        eMatK = ElementalKinectMatrix(order);

        for(int ei=0; ei<Ne; ei++)
        {
                coeff = 0.5*e[ei].h;
		vij[0] = e[ei].pot[0].x;
		vij[1] = e[ei].pot[1].x;
		vij[2] = e[ei].pot[2].x;
		vij[3] = e[ei].pot[3].x;
		eMatV = ElementalPotentialMatrix(order,vij);
                for(int nu=0; nu<p; nu++)
                {
                        l=link_mat[p*ei+nu];
                        for(int mu=0; mu<p; mu++)
                        {
                                m=link_mat[p*ei+mu];
				s_mat[l*nodes+m] +=  coeff*eMatS[p*nu+mu];
				k_mat[l*nodes+m] +=  eMatK[p*nu+mu]/coeff;
				v_mat[l*nodes+m] +=  coeff*eMatV[p*nu+mu];
                                //printf("%d %d\n",l*nodos+m,p*nu+mu);
                        }
                }
        }
	free(eMatV);
	free(eMatK);
	free(eMatS);


}
void AssambleGlobalMatrices(int Ne,int order,int *link_mat,double *sij,double *kij,double *vij,double *uij,double *lij,struct Element *e,struct Vertex *pot)
{
	int p,l,m;
	int nodes=(order*Ne) + 1;
	p = order+1;
	double *s_mat = (double *)malloc(sizeof(double)*(nodes*nodes));
	double *k_mat = (double *)malloc(sizeof(double)*(nodes*nodes));
	double *v_mat = (double *)malloc(sizeof(double)*(nodes*nodes));
	double *u_mat = (double *)malloc(sizeof(double)*(nodes*nodes));
	double *l_mat = (double *)malloc(sizeof(double)*(nodes*nodes));
	FillZeroMat(s_mat,nodes,nodes);
	FillZeroMat(k_mat,nodes,nodes);
	FillZeroMat(v_mat,nodes,nodes);
	FillZeroMat(u_mat,nodes,nodes);
	FillZeroMat(l_mat,nodes,nodes);
	switch(order)
	{
		case 1:
			FirstGradePol(Ne,order,link_mat,s_mat,k_mat,v_mat,e,e->pot);
			break;
		case 2: 
			SecondGradePol(Ne,order,link_mat,s_mat,k_mat,v_mat,e,e->pot);
			break;
		case 3:
			ThirdGradePol(Ne,order,link_mat,s_mat,k_mat,v_mat,e,e->pot);
			break;
	}

	ScalarXMatrix(1.0,s_mat,u_mat,nodes,nodes);
	ScalarXMatrix(2.0,k_mat,l_mat,nodes,nodes);


	ApplyBoundaryCond(s_mat,k_mat,v_mat,u_mat,l_mat,nodes);
	ReduceMatrices(sij,kij,vij,uij,lij,s_mat,k_mat,v_mat,u_mat,l_mat,nodes);
	/*print_matrix("*** Global Overlap Matrix ***",nodes,nodes,s_mat);
	print_matrix("*** Global Kinect Matrix ***",nodes,nodes,k_mat);
        print_matrix("*** Global Potential Matrix ***",nodes,nodes,v_mat);
        print_matrix("*** Global U Poisson Matrix ***",nodes,nodes,u_mat);
        print_matrix("*** Global L Poisson Matrix ***",nodes,nodes,l_mat);*/

	free(l_mat);
	free(u_mat);
	free(v_mat);
	free(k_mat);
	free(s_mat);


}
