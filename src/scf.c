#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "data_structures.h"
	
extern void FillZeroMat(double *mat,int M,int N);
extern void NormWfn(double *wfn_vec,int *link_mat,struct Element *e,int order,int Ne);
extern void GetHartreePotential(double *mat_vh,double *sij,double *kij,double *wfn,int Ne,int order);
extern void diag (int n, double *h, double *s, double *e, double *v);
extern void GetHmatrix(double *h,double *v,double *k,int M,int N);

void PerformSCF(double *hij,double *sij,double *kij,double *vij,double *vh_mat,double *wfn,double *vec_ei,int *link_mat,struct Element *e,int Ne,int order)
{
        printf("****************************************************************************\n");
	printf("**************************  SCF MODULE  ************************************\n");
        printf("****************************************************************************\n");
        printf("****************************************************************************\n");
	int nodes = Ne*order + 1;
	int r_nodes = nodes-2;
	int K=r_nodes;

	
	double *f_mat = (double *)malloc(sizeof(double)*(K*K));

	//   H_ij = K_ij + VN_ij + VEE_ij
	//FIRST APPROACH TO H_ij MATRIX: 
	//   Hij = K_ij + VN_ij
	FillZeroMat(vh_mat,r_nodes,r_nodes);
	FillZeroMat(hij,r_nodes,r_nodes);
	GetHmatrix(hij,vij,kij,r_nodes,r_nodes);
	int count=0;
	double energy = 0.0;
	double orbE_new,orbE_old,delta_energy;
	delta_energy=10000000000.0;
	orbE_old=0.0;
	do
	{
		printf("\n");
		for(int i=0; i<r_nodes*r_nodes; i++)
		{
			f_mat[i] = hij[i] + vh_mat[i];
		}	
		diag(K,f_mat,sij,vec_ei,wfn);
		orbE_new = 0.5*vec_ei[0];
		printf("Orbital[%d] energy = %lf    Step: %d\n",0,orbE_new,count);
		NormWfn(wfn,link_mat,e,order,Ne);
		GetHartreePotential(vh_mat,sij,kij,wfn,Ne,order);
		delta_energy = fabs(orbE_old-orbE_new);
		orbE_old = orbE_new;
		/*for(int i=0; i<r_nodes*r_nodes; i++)
		{
			energy = energy + 0.5*hij[i]*wfn[i];
			energy = energy + 0.5*f_mat[i]*wfn[i];
		}*/

		//printf("TOTAL ENERGY = %lf\n",energy);



		count++;
	}while(delta_energy>0.00000001);
        printf("****************************************************************************\n");

}
