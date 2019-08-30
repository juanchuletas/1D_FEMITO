#include <stdio.h>
#include <stdlib.h>
#include "wfn_operations.h"
extern void dgesv_( int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info );
void GetHartreePotential(double *mat_vh,double *sij,double *kij,double *wfn,int Ne,int order,int phase)
{
	int nodes = Ne*order + 1;
	int r_nodes = nodes-2;
	int N=r_nodes;
	int NRHS=r_nodes;
	int LDA=N;
	int LDB=N;
	int ipiv[N];
	int info;
	int orb=0;
	double rho,coeffwfn;
	const double PI = 3.14159265358979323846;

	//LOCAL ARRAYS:
	//double vec_h[r_nodes*r_nodes],aux_mat[r_nodes*r_nodes],mat_kij[r_nodes*r_nodes];
	double vec_h[r_nodes],aux_mat[r_nodes*r_nodes],mat_kij[r_nodes*r_nodes];
	int k=0;
	for(int i=0; i<r_nodes; i++)
	{
		coeffwfn = wfn[i + orb*r_nodes]*phase;
		rho = coeffwfn*coeffwfn;
		for(int j=0; j<r_nodes; j++)
		{
			aux_mat[k] = sij[k]*(rho);
			mat_kij[k] = 2.0*kij[k];
			k++;

		}
	}
	dgesv_(&N,&NRHS,mat_kij,&LDA,ipiv,aux_mat,&LDB,&info);
	if( info > 0 ) 
	{
        printf( "The diagonal element of the triangular factor of A,\n" );
        printf( "U(%i,%i) is zero, so that A is singular;\n", info, info );
        printf( "the solution could not be computed.\n" );
        exit(1);
	}

	for(int i=0; i<r_nodes*r_nodes; i++)
	{
		mat_vh[i] = aux_mat[i];
		//printf("Hartree Pot = %lf\n",mat_vh[i]);
	}

}
