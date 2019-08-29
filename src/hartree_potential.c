#include <stdio.h>
#include <stdlib.h>
extern void dgesv_( int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info );
void GetHartreePotential(double *mat_vh,double *sij,double *kij,double *wfn,int Ne,int order)
{
	int nodes = Ne*order + 1;
	int r_nodes = nodes-2;
	int N=r_nodes;
	int NRHS=r_nodes;
	int LDA=N;
	int LDB=N;
	int ipiv[N];
	int info;
	int orb = r_nodes;
	//LOCAL ARRAYS:
	double vec_h[r_nodes*r_nodes],aux_mat[r_nodes*r_nodes],mat_kij[r_nodes*r_nodes];

	for(int i=0; i<r_nodes*r_nodes; i++)
	{
		aux_mat[i] = sij[i]*(wfn[i]*wfn[i]);
		mat_kij[i] = 2.0*kij[i];
		//printf("2kij=%lf   Sij=%lf  rho=%lf\n",mat_kij[i],sij[i],wfn[i]*wfn[i]);
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
		
	//	printf("Hartree Pot = %lf\n",mat_vh[i]);
	}

}
