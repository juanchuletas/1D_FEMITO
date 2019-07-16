#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "data_structures.h"


extern void GenerateNodes(double *x,int Ne,double r0,double rN,char mesh[180]);
extern void SetElements(struct Element *e,struct Vertex *n,int Ne, double *x,int order);
extern void GetElementalMatrices(int order,double *eMatS,double *eMatK,double *eMatV,double *v);
extern void EvaluatePotential(char pot[180],int nodes,double *x,double *v);
extern void print_matrix( char* desc, int m, int n, double* a);

void LinearInterpolation(int Ne,double r0,double rN,char pot[180],char mesh[180],int order)
{

	int M,me,N,ne,nodes;
	nodes = Ne+1;
	me = (order+1);
	me = (order+1);
	M = (Ne+1);
	N = (Ne+1);
	int eMatSize = ne*ne;
	double *s_mat,*k_mat,*h_mat,*v_mat,*eMatS;
	double *eMatK,*eMatV,*x,*se_mat,*v;

	struct Element *e;

	//Structures for nodes********************************
        e = (struct Element *)malloc(sizeof(struct Element)*Ne);
	//Points in x ****************************************
	x = (double *)malloc(sizeof(double)*(nodes));
	v = (double *)malloc(sizeof(double)*(nodes));
	// Memory allocation for the elemental matrices******
	eMatS = (double *)malloc(sizeof(double)*(eMatSize));
	eMatK = (double *)malloc(sizeof(double)*(eMatSize));
	eMatV = (double *)malloc(sizeof(double)*(eMatSize));
	// Memory allocation for the global matrices*********
	
	/*s_mat = (double *)malloc(sizeof(double)*((Ne+1)*(Ne+1)));
	k_mat = (double *)malloc(sizeof(double)*((Ne+1)*(Ne+1)));
	v_mat = (double *)malloc(sizeof(double)*((Ne+1)*(Ne+1)));
	h_mat = (double *)malloc(sizeof(double)*((Ne+1)*(Ne+1)));*/
	// Generation for the nodes and the elements ***********
	GenerateNodes(x,Ne,r0,rN,mesh);
	SetElements(e,e->n,Ne,x,order);
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

	/*SetZeroMat(s_mat,M,N);	// FILL WITH 0.0 THE MATRIX S
	SetZeroMat(k_mat,M,N);	// FILL WITH 0.0 THE MATRIX K
	SetZeroMat(v_mat,M,N);	// FILL WITH 0.0 THE MATRIX V
	SetZeroMat(h_mat,M,N);	// FILL WITH 0.0 THE MATRIX H*/

	
	EvaluatePotential(pot,nodes,x,v);

	GetElementalMatrices(order,eMatS,eMatK,eMatV,v);
	print_matrix("Elemental Overlap Matrix",me,ne,eMatS);
	print_matrix("Elemental Kinect Matrix",me,ne,eMatK);
	print_matrix("Elemental Potential Matrix",me,ne,eMatV);


	free(eMatS);
	free(eMatK);
	free(eMatV);

}
