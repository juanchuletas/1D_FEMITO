#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "data_structures.h"
#include "wfn_operations.h"
extern void FillZeroMat(double *mat,int M,int N);
extern double *PoissonSolver(double *uij,double *lij,double *wfn,int Ne,int order);
extern void AssambleHartreePot(int Ne,int order,int *link_mat,double *vhij,struct Element *e,double *v);
extern void print_matrix(char *name,int m,int n,double *matA);
extern void NormWfn(double *wfn_vec,int *link_mat,struct Element *e,int order,int Ne);
extern void GetHartreePotential(double *mat_vh,double *sij,double *kij,double *wfn,double *vec_ei,int Ne,int order,int phase);
extern void diag (int n, double *h, double *s, double *e, double *v);
extern void GetHmatrix(double *h,double *v,double *k,int M,int N);
extern void ComputeDensityMatrix(double *dens_mat,double *wfn,int Ne, int order, int atomicN);
void PerformSCF(double *hij,double *sij,double *kij,double *vij,double *uij,double *lij,double *wfn,double *vec_ei,int *link_mat,struct Element *e,int Ne,int order)
{
        printf("****************************************************************************\n");
	printf("**************************  SCF MODULE  ************************************\n");
        printf("****************************************************************************\n");
        printf("****************************************************************************\n");
	int nodes = Ne*order + 1;
	int r_nodes = nodes-2;
	int K=r_nodes;
	int orb=0;
	
	double *f_mat = (double *)malloc(sizeof(double)*(K*K));
	double *hartree_vec = (double *)malloc(sizeof(double)*(nodes));
	double *dens_mat = (double *)malloc(sizeof(double)*(K*K));
	double *vhij = (double *)malloc(sizeof(double)*(K*K));

	//   H_ij = K_ij + VN_ij + VEE_ij
	//FIRST APPROACH TO H_ij MATRIX: 
	//   Hij = K_ij + VN_ij
	FillZeroMat(hij,r_nodes,r_nodes);
	GetHmatrix(hij,vij,kij,r_nodes,r_nodes);
	int count=0;
	int phase=1;
	double energy = 0.0;
	double orbE_new,orbE_old,delta_energy;
	delta_energy=10000000000.0;
	orbE_old=0.0;
	//do
	//{
		printf("\n");
		/*for(int i=0; i<r_nodes*r_nodes; i++)
		{
			f_mat[i] = hij[i];
		}*/	
		diag(K,hij,sij,vec_ei,wfn);
		orbE_new = vec_ei[0];
		printf("Orbital[%d] energy = %lf    Step: %d\n",0,orbE_new,count);
		NormWfn(wfn,link_mat,e,order,Ne);
		GetWfnPhase(r_nodes,orb,&phase,wfn);
		hartree_vec = PoissonSolver(uij,lij,wfn,Ne,order);
		printf("HARTREE POTENTIAL AT NODES:\n");
		for(int i=0; i<nodes; i++)
		{
			printf("%d   %lf\n",i,hartree_vec[i]);
		}
		AssambleHartreePot(Ne,order,link_mat,vhij,e,hartree_vec);
		print_matrix("*** Hartree Matrix ***",K,K,vhij);

		//ComputeDensityMatrix(dens_mat,wfn,Ne,order,2);

		//count++;
	//}while(count<30);
        printf("****************************************************************************\n");

	free(f_mat);
	free(dens_mat);
	free(vhij);

}
