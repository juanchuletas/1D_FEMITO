#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "data_structures.h"
#include <math.h>

extern void GenerateNodes(double *x,int nodes,double r0,double rN,char mesh[180]);
extern void SetElements(struct Element *e,struct Vertex *n,int Ne, double *x,int order,double *nt);
extern void EvaluatePotential(char kindpot[180],int Ne,int order,struct Element *e,struct Vertex *n,struct Vertex *pot);
extern void print_matrix( char* desc, int m, int n, double* a);
extern void GetLinkMatrix(int *link_mat,int Ne,int order);
extern void AssambleGlobalMatrices(int Ne,int order,int *link_mat,double *s_mat,double *k_mat,double *v_mat,double *v,struct Element *e,struct Vertex *pot);
extern void FillZeroMat(double *mat,int M,int N);
extern void GetHmatrix(double *h,double *v,double *k,int M,int N);
extern void diag (int n, double *h, double *s, double *e, double *v);
extern void ReduceMatrices(double *sij,double *kij,double *vij,double *s_mat,double *k_mat,double *v_mat,int nodes);

void LagrangeInterpolation(int Ne,double r0,double rN,char kindpot[180],char mesh[180],int order)
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
	double *x,*v,*nt;
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
	nt = (double *)malloc(sizeof(double)*(nodes));
	v = (double *)malloc(sizeof(double)*(nodes));
	//***************** Link Matrix ****************************************
	link_mat = (int *)malloc(sizeof(int)*(Ne*(order+1)));
	//***************** Eigenvectos and eigenvalues ************************
	ci = (double *)malloc(pow(M,2)*sizeof(double));
        ei = (double *)malloc((N)*sizeof(double));

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
	ci = (double *)malloc(pow(K,2)*sizeof(double));
        ei = (double *)malloc((K)*sizeof(double));
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------



	//------------ GENERATION OF NODES AND ELEMENTS------------------------
	GenerateNodes(x,Ne,r0,rN,mesh);
	SetElements(e,e->n,Ne,x,order,nt);
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
	EvaluatePotential(kindpot,Ne,order,e,e->n,e->pot);
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
	print_matrix("EIGENVALUES",1,K,ei);
	//Perform_SCF(sij,kij,vij,nodes);

	/*for(int i=0; i<Ne; i++)
        {
                        printf("e[%d].h = %lf\n",i,e[i].h);
        }
	printf("NODE PER ELEMENT VALUE:\n");
        for(int i=0; i<Ne; i++)
        {
                for(int j=0; j<p; j++)
                {
                        printf("e[%d].n[%d].x = %lf\n",i,j,e[i].n[j].x);

                }
        }
	printf("NODE PER ELEMENT VALUE:\n");
        for(int i=0; i<Ne; i++)
        {
                for(int j=0; j<p; j++)
                {
                        printf("e[%d].v[%d].x = %lf\n",i,j,e[i].pot[j].x);

                }
        }*/


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
