#include <stdio.h>
#include <stdlib.h>
extern void MatrixProduct(double *A,double *B,double *C,int N, int M, int P);
extern void ScalarXMatrix(double coeff,double *mat_A,double *mat_Res,int N,int M);
extern void dgesv_( int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info );
void GetHartreePotential(double *hartree_vec,double *sij,double *kij,double *wfn,int Ne,int order,int phase)
{
	int nodes = Ne*order + 1;
	int r_nodes = nodes-2;
	int N=r_nodes;
	int NRHS=1;
	int LDA=N;
	int LDB=N;
	int ipiv[N];
	int info;
	int orb=0;
	double rho[r_nodes],coeffwfn;
	const double PI = 3.14159265358979323846;

	//LOCAL ARRAYS:
	//double vec_h[r_nodes*r_nodes],aux_mat[r_nodes*r_nodes],mat_kij[r_nodes*r_nodes];
	double aux_mat[r_nodes*r_nodes],aux_kij[r_nodes*r_nodes];
	double right_vec[r_nodes],aux_vec[r_nodes];

	for(int i=0; i<r_nodes; i++)
	{
		coeffwfn = wfn[i + orb*r_nodes]*phase;
		rho[i] = coeffwfn*coeffwfn;
	}

	MatrixProduct(sij,rho,right_vec,r_nodes,r_nodes,NRHS);
	ScalarXMatrix(2.0,kij,aux_kij,r_nodes,r_nodes);
	ScalarXMatrix(4.0*PI,right_vec,aux_vec,r_nodes,NRHS);


	dgesv_(&N,&NRHS,aux_kij,&LDA,ipiv,aux_vec,&LDB,&info);

	if( info > 0 ) 
	{
        printf( "The diagonal element of the triangular factor of A,\n" );
        printf( "U(%i,%i) is zero, so that A is singular;\n", info, info );
        printf( "the solution could not be computed.\n" );
        exit(1);
	}
	for(int i=0; i<r_nodes; i++)
	{
		hartree_vec[i] = aux_vec[i];
		printf("Hartree Potiential = %lf\n",hartree_vec[i]);
	}

}
