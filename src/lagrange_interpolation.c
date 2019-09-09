#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "data_structures.h"
#include <math.h>
extern void PrintWaveFunction(int Ne,int order,struct Element *e,struct Vertex *n,double *wfn,double *vecE);
extern void EnergyResults(double *eigenval,int N);
extern void PrintPotentialData(int order,int Ne,struct Element *e,struct Vertex *n,struct Vertex *pot);
extern void GenerateNodes(double *x,int nodes,double r0,double rN,char mesh[180],char atom[3]);
extern void SetElements(struct Element *e,struct Vertex *n,int Ne, double *x,int order);
extern void EvaluatePotential(char kindpot[180],char atom[3],int angular,int Ne,int order,struct Element *e,struct Vertex *n,struct Vertex *pot);
extern void print_matrix( char* desc, int m, int n, double* a);
extern void GetLinkMatrix(int *link_mat,int Ne,int order);
extern void AssambleGlobalMatrices(int Ne,int order,int *link_mat,double *sij,double *kij,double *vij,double *uij,double *lij,struct Element *e,struct Vertex *pot);
//extern void FillZeroMat(double *mat,int M,int N);
//extern void GetHmatrix(double *h,double *v,double *k,int M,int N);
//extern void diag (int n, double *h, double *s, double *e, double *v);
//extern void ReduceMatrices(double *sij,double *kij,double *vij,double *s_mat,double *k_mat,double *v_mat,int nodes);
//extern void NormWfn(double *wfn_vec,int *link_mat,struct Element *e,int order,int Ne);
//extern void GetHartreePotential(double *mat_vh,double *sij,double *kij,double *wfn,int Ne,int order);
extern void PerformSCF(double *hij,double *sij,double *kij,double *vij,double *uij,double *lij,double *wfn,double *vec_ei,int *link_mat,struct Element *e,int Ne,int order);

void LagrangeInterpolation(int Ne,double r0,double rN,char kindpot[180],char mesh[180],char atom[3],int order,int angular)
{
	int nodes = Ne*order + 1;
	int fembasis = nodes-2;
	int fembasis_poisson = nodes-1;

	int *link_mat;
	struct Element *e;
	double *x;

	//---------- MEMORY ALLOCATION FOR THE NECCESARY MATRICES ----------
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	//******************Structures for nodes********************************
        e = (struct Element *)malloc(sizeof(struct Element)*Ne);
	//***************** Nodes in x ****************************************
	x = (double *)malloc(sizeof(double)*(Ne+1));
	//***************** Link Matrix ****************************************
	link_mat = (int *)malloc(sizeof(int)*(Ne*(order+1)));
	
	// **** Memory allocation: Reduced matrices for eigenvalue problem
	double *sij,*kij,*vij,*hij,*vh_mat,*uij,*lij,*ci,*ei,*vhij;
	int K = nodes-2;
	sij = (double *)malloc(sizeof(double)*(fembasis*fembasis));
	kij = (double *)malloc(sizeof(double)*(fembasis*fembasis));
	vij = (double *)malloc(sizeof(double)*(fembasis*fembasis));
	hij = (double *)malloc(sizeof(double)*(fembasis*fembasis));
	vhij = (double *)malloc(sizeof(double)*(fembasis*fembasis));
	//***************** Eigenvectos and eigenvalues ************************
	ci = (double *)malloc(sizeof(double)*(K*K));
        ei = (double *)malloc((K)*sizeof(double));
	// **** Memory allocation: Reduced matrices for the Poisson problem
	lij = (double *)malloc(sizeof(double)*(fembasis_poisson*fembasis_poisson));
	uij = (double *)malloc(sizeof(double)*(fembasis_poisson*fembasis_poisson));
	

	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------



	//------------ GENERATION OF NODES AND ELEMENTS------------------------
	GenerateNodes(x,Ne,r0,rN,mesh,atom);
	SetElements(e,e->n,Ne,x,order);
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	//------------------- EVALUATE THE POTENTIAL ---------------------------
	//
	EvaluatePotential(kindpot,atom,angular,Ne,order,e,e->n,e->pot);
	//
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------

	//------------------- GET THE ELEMENTAL MATRICES AND -------------------
	//--------------------------- THE LINK MATRIX ---------------------------
	GetLinkMatrix(link_mat,Ne,order);
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	//------------------- GET THE GLOBAL MATRICES ---------------------------
	AssambleGlobalMatrices(Ne,order,link_mat,sij,kij,vij,uij,lij,e,e->pot);
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	print_matrix("*** Overlap Matrix ***",K,K,sij);
	print_matrix("*** Kinect Matrix ***",K,K,kij);
        print_matrix("*** Potential Matrix ***",K,K,vij);
        print_matrix("*** U Poisson Matrix ***",fembasis_poisson,fembasis_poisson,uij);
        print_matrix("*** L Poisson Matrix ***",fembasis_poisson,fembasis_poisson,lij);
	PrintPotentialData(order,Ne,e,e->n,e->pot);
        PerformSCF(hij,sij,kij,vij,uij,lij,ci,ei,link_mat,e,Ne,order);
	//EnergyResults(ei,10);
	//print_matrix("Hartree-Potential",K,K,mat_vh);
	PrintWaveFunction(Ne,order,e,e->n,ci,ei);

	free(link_mat);
	free(hij);
	free(vhij);
	free(sij);
	free(kij);
	free(uij);
	free(lij);
	free(vij);
	free(e);
	free(x);
	free(ei);
	free(ci);

}
