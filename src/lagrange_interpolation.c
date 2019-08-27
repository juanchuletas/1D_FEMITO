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
extern void AssambleGlobalMatrices(int Ne,int order,int *link_mat,double *s_mat,double *k_mat,double *v_mat,double *v,struct Element *e,struct Vertex *pot);
extern void FillZeroMat(double *mat,int M,int N);
extern void GetHmatrix(double *h,double *v,double *k,int M,int N);
extern void diag (int n, double *h, double *s, double *e, double *v);
extern void ReduceMatrices(double *sij,double *kij,double *vij,double *s_mat,double *k_mat,double *v_mat,int nodes);
extern void NormWfn(double *wfn_vec,int *link_mat,struct Element *e,int order,int Ne);
extern int normWF(int nele,int p,int *g, struct Element *e, double *matC);
void LagrangeInterpolation(int Ne,double r0,double rN,char kindpot[180],char mesh[180],char atom[3],int order,int angular)
{

	int M,me,N,ne,nodes,p;
	p = order + 1;
	nodes = (order*Ne) + 1;
	me = (order+1);
	ne = (order+1);
	M = (Ne+1);
	N = (Ne+1);
	int eMatSize = ne*me;
	double *s_mat,*k_mat,*v_mat,*ei,*ci;
	double *x,*v;
	int *link_mat;
	struct Element *e;

	//---------- MEMORY ALLOCATION FOR THE NECCESARY MATRICES ----------
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	//******************Structures for nodes********************************
        e = (struct Element *)malloc(sizeof(struct Element)*Ne);
	//***************** Nodes in x ****************************************
	x = (double *)malloc(sizeof(double)*(Ne+1));
	//***************** Potential ******************************************
	v = (double *)malloc(sizeof(double)*(nodes));
	//***************** Link Matrix ****************************************
	link_mat = (int *)malloc(sizeof(int)*(Ne*(order+1)));
	// Memory allocation for the elemental matrices******
	//eMatS = (double *)malloc(sizeof(double)*(eMatSize));
	//eMatK = (double *)malloc(sizeof(double)*(eMatSize));
	//eMatV = (double *)malloc(sizeof(double)*(eMatSize));
	// Memory allocation for the global matrices*********
	
	s_mat = (double *)malloc(sizeof(double)*(nodes*nodes));
	k_mat = (double *)malloc(sizeof(double)*(nodes*nodes));
	v_mat = (double *)malloc(sizeof(double)*(nodes*nodes));
	//
	// **** Memory allocation for reduces matrices
	double *sij,*kij,*vij,*hij;
	int K = nodes-2;
	sij = (double *)malloc(sizeof(double)*(K*K));
	kij = (double *)malloc(sizeof(double)*(K*K));
	vij = (double *)malloc(sizeof(double)*(K*K));
	hij = (double *)malloc(sizeof(double)*(K*K));
	//***************** Eigenvectos and eigenvalues ************************
	ci = (double *)malloc(sizeof(double)*(K*K));
        ei = (double *)malloc((K)*sizeof(double));
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
	//------------- FILL THE GLOBAL MATRICES WITH ZEROS --------------------
	FillZeroMat(s_mat,nodes,nodes);	// FILL WITH 0.0 THE MATRIX S
	FillZeroMat(k_mat,nodes,nodes);	// FILL WITH 0.0 THE MATRIX K
	FillZeroMat(v_mat,nodes,nodes);	// FILL WITH 0.0 THE MATRIX V
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
	AssambleGlobalMatrices(Ne,order,link_mat,s_mat,k_mat,v_mat,v,e,e->pot);
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	//print_matrix("*** Overlap Matrix ***",M,N,s_mat);
	//print_matrix("*** Kinect Matrix ***",M,N,k_mat);
        //print_matrix("*** Potential Matrix ***",M,N,v_mat);
	//print_matrix("*** H Matrix ***",M,N,h_mat);
	ReduceMatrices(sij,kij,vij,s_mat,k_mat,v_mat,nodes);
	//print_matrix("*** Overlap Matrix ***",K,K,sij);
	//print_matrix("*** Kinect Matrix ***",K,K,kij);
	//print_matrix("*** Potential Matrix ***",K,K,vij);
	//---------------------------------------------------
	FillZeroMat(hij,K,K);	// FILL WITH 0.0 THE MATRIX V
	GetHmatrix(hij,vij,kij,K,K);
	//print_matrix("*** H  Matrix ***",K,K,hij);
	diag(K,hij,sij,ei,ci);
	//print_matrix("EIGENVALUES",1,K,ei);
	PrintPotentialData(order,Ne,e,e->n,e->pot);
	EnergyResults(ei,10);
	PrintPotentialData(order,Ne,e,e->n,e->pot);
	//normWF(Ne,order,link_mat,e,ci);
	NormWfn(ci,link_mat,e,order,Ne);
	//print_matrix("EIGENVECTORS",K,K,ci);
	PrintWaveFunction(Ne,order,e,e->n,ci,ei);
	//Perform_SCF(sij,kij,vij,nodes);

	free(k_mat);
	free(v_mat);
	free(s_mat);
	free(link_mat);
	free(hij);
	free(sij);
	free(kij);
	free(vij);
	free(e);
	free(x);
	free(ei);
	free(ci);

}
