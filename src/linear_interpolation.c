#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "data_structures.h"
#include <math.h>

extern void GenerateNodes(double *x,int Ne,double r0,double rN,char mesh[180]);
extern void SetElements(struct Element *e,struct Vertex *n,int Ne, double *x,int order);
extern void GetElementalMatrices(int order,double *eMatS,double *eMatK,double *eMatV,double *v);
extern void EvaluatePotential(char pot[180],int nodes,double *x,double *v);
extern void print_matrix( char* desc, int m, int n, double* a);
extern void GetLinkMatrix(double *link_mat,int Ne,int order);
extern void AssambleGlobalMatrices(int Ne,int order,double *link_mat,double *s_mat,double *k_mat,double *v_mat,double *eMatS,double *eMatK,double *eMatV,struct Element *e);
extern void FillZeroMat(double *mat,int M,int N);
extern void GetHmatrix(double *h,double *v,double *k,int M,int N);
extern void diag (int n, double *h, double *s, double *e, double *v);
extern void ReduceMatrices(double *sij,double *kij,double *vij,double *s_mat,double *k_mat,double *v_mat,int nodes);

void LinearInterpolation(int Ne,double r0,double rN,char pot[180],char mesh[180],int order)
{

	int M,me,N,ne,nodes;
	nodes = Ne+1;
	me = (order+1);
	ne = (order+1);
	M = (Ne+1);
	N = (Ne+1);
	int eMatSize = ne*me;
	double *s_mat,*k_mat,*v_mat,*eMatS,*ei,*ci;
	double *eMatK,*eMatV,*x,*v,*link_mat;
	struct Element *e;

	//---------- MEMORY ALLOCATION FOR THE NECCESARY MATRICES ----------
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	//******************Structures for nodes********************************
        e = (struct Element *)malloc(sizeof(struct Element)*Ne);
	//***************** Nodes in x ****************************************
	x = (double *)malloc(sizeof(double)*(nodes));
	//***************** Potential ******************************************
	v = (double *)malloc(sizeof(double)*(nodes));
	//***************** Link Matrix ****************************************
	link_mat = (double *)malloc(sizeof(double)*(Ne*(order+1)));
	//***************** Eigenvectos and eigenvalues ************************
	ci = (double *)malloc(pow(M,2)*sizeof(double));
        ei = (double *)malloc((N)*sizeof(double));

	// Memory allocation for the elemental matrices******
	eMatS = (double *)malloc(sizeof(double)*(eMatSize));
	eMatK = (double *)malloc(sizeof(double)*(eMatSize));
	eMatV = (double *)malloc(sizeof(double)*(eMatSize));
	// Memory allocation for the global matrices*********
	
	s_mat = (double *)malloc(sizeof(double)*((Ne+1)*(Ne+1)));
	k_mat = (double *)malloc(sizeof(double)*((Ne+1)*(Ne+1)));
	v_mat = (double *)malloc(sizeof(double)*((Ne+1)*(Ne+1)));
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
	SetElements(e,e->n,Ne,x,order);
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	for(int i=0; i<Ne; i++)
        {
                        printf("e[%d].h = %lf\n",i,e[i].h);
        }
        for(int i=0; i<Ne; i++)
        {
                for(int j=0; j<2; j++)
                {
                        printf("e[%d].n[%d].x = %lf\n",i,j,e[i].n[j].x);

                }
        }

	//------------- FILL THE GLOBAL MATRICES WITH ZEROS --------------------
	FillZeroMat(s_mat,M,N);	// FILL WITH 0.0 THE MATRIX S
	FillZeroMat(k_mat,M,N);	// FILL WITH 0.0 THE MATRIX K
	FillZeroMat(v_mat,M,N);	// FILL WITH 0.0 THE MATRIX V
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	//------------------- EVALUATE THE POTENTIAL ---------------------------
	//
	EvaluatePotential(pot,nodes,x,v);
	//
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------

	//------------------- GET THE ELEMENTAL MATRICES AND -------------------
	//--------------------------- THE LINK MATRIX ---------------------------
	GetElementalMatrices(order,eMatS,eMatK,eMatV,v);
	GetLinkMatrix(link_mat,Ne,order);
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	print_matrix("*** Elemental Overlap Matrix ***",me,ne,eMatS);
	print_matrix("*** Elemental Kinect Matrix ***",me,ne,eMatK);
	print_matrix("*** Elemental Potential Matrix ***",me,ne,eMatV);
	//------------------- GET THE GLOBAL MATRICES ---------------------------
	AssambleGlobalMatrices(Ne,order,link_mat,s_mat,k_mat,v_mat,eMatS,eMatK,eMatV,e);
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	//print_matrix("*** Overlap Matrix ***",M,N,s_mat);
	//print_matrix("*** Kinect Matrix ***",M,N,k_mat);
        //print_matrix("*** Potential Matrix ***",M,N,v_mat);
	//print_matrix("*** H Matrix ***",M,N,h_mat);
	ReduceMatrices(sij,kij,vij,s_mat,k_mat,v_mat,nodes);
	print_matrix("*** Overlap Matrix ***",K,K,sij);
	print_matrix("*** Kinect Matrix ***",K,K,kij);
	print_matrix("*** Potential Matrix ***",K,K,vij);
	//---------------------------------------------------
	//FillZeroMat(hij,K,K);	// FILL WITH 0.0 THE MATRIX V
	//GetHmatrix(hij,vij,kij,K,K);
	//print_matrix("*** H  Matrix ***",K,K,hij);
	//diag(K,hij,sij,ei,ci);
	//print_matrix("EIGENVALUES",1,K,ei);
	//Perform_SCF(sij,kij,vij,nodes);


	free(eMatS);
	free(eMatK);
	free(eMatV);
	free(k_mat);
	free(v_mat);
	free(s_mat);
	free(hij);
	free(sij);
	free(kij);
	free(vij);
	free(e);
	free(x);
	free(ei);
	free(ci);

}
