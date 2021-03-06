#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "data_structures.h"
#include "wfn_operations.h"
extern double *ComputeDensityR(double *wfn,struct Element *e,struct Vertex *n,int Ne, int order, int atomicN,int sign);
extern void FillZeroMat(double *mat,int M,int N);
extern double *PoissonSolver(double *lij,double *wfn,int *link_mat,struct Element *e,int Ne,int order,int phase);
extern void AssambleHartreePot(int Ne,int order,int *link_mat,double *vhij,struct Element *e,double *v);
extern void print_matrix(char *name,int m,int n,double *matA);
extern void NormWfn(double *wfn_vec,int *link_mat,struct Element *e,int order,int Ne);
extern void GetHartreePotential(double *mat_vh,double *sij,double *kij,double *wfn,double *vec_ei,int Ne,int order,int phase);
extern void diag (int n, double *h, double *s, double *e, double *v);
extern void GetHmatrix(double *h,double *v,double *k,int M,int N);
extern double *ComputeDensity(double *wfn,int Ne, int order, int atomicN);
extern double IntegrateDens(double *wfn,struct Element *e,int *link_mat,int Ne, int order);
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
	double *hartree_vec = (double *)malloc(sizeof(double)*(r_nodes));
	double *rho = (double *)malloc(sizeof(double)*(r_nodes));
	double *vhij = (double *)malloc(sizeof(double)*(K*K));

	//   H_ij = K_ij + VN_ij + VEE_ij
	//FIRST APPROACH TO H_ij MATRIX: 
	//   Hij = K_ij + VN_ij
	FillZeroMat(hij,r_nodes,r_nodes);
	FillZeroMat(f_mat,r_nodes,r_nodes);
	GetHmatrix(hij,vij,kij,r_nodes,r_nodes);
	int count=0;
	int phase=1;
	double energy = 0.0;
	double orbE_new,orbE_old,delta_energy;
	delta_energy=10000000000.0;
	orbE_old=0.0;



	diag(K,hij,sij,vec_ei,wfn);
        orbE_new = vec_ei[0];
        printf("Orbital[%d] energy = %lf    Step: %d\n",0,orbE_new,count);
        NormWfn(wfn,link_mat,e,order,Ne);
        GetWfnPhase(r_nodes,orb,&phase,wfn);
	double qtot = IntegrateDens(wfn,e,link_mat,Ne,order);

	printf("TOTAL CHARGE = %lf\n",qtot);

	rho = ComputeDensity(wfn,Ne,order,2);
	/*for(int i=0; i<r_nodes; i++)
	{
		printf("RHO[%d] = %lf\n",i,rho[i]);
	}*/
        //hartree_vec = PoissonSolver(kij,rho,link_mat,e,Ne,order,phase);
	//AssambleHartreePot(Ne,order,link_mat,vhij,e,hartree_vec);
	/*do
	{
	        GetHmatrix(f_mat,hij,vhij,r_nodes,r_nodes);
		diag(K,f_mat,sij,vec_ei,wfn);
        	printf("Orbital[%d] energy = %lf    Step: %d\n",0,vec_ei[0],count+1);
        	//NormWfn(wfn,link_mat,e,order,Ne);
        	//GetWfnPhase(r_nodes,orb,&phase,wfn);
	        rho = ComputeDensity(wfn,Ne,order,2);
        	hartree_vec = PoissonSolver(kij,rho,link_mat,e,Ne,order,phase);
		AssambleHartreePot(Ne,order,link_mat,vhij,e,hartree_vec);
		count++;

	}while(count<100);
	//ComputeDensityMatrix(dens_mat,wfn,Ne,order,2);

        printf("****************************************************************************\n");
*/
	free(f_mat);
	free(rho);
	free(vhij);
	free(hartree_vec);

}
