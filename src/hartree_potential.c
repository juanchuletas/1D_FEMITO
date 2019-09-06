#include <stdio.h>
#include <stdlib.h>
extern void ColumnMayor(double *mat_A,double *mat_C,int N, int M);
extern void MatrixProduct(double *A,double *B,double *C,int N, int M, int P);
extern void ScalarXMatrix(double coeff,double *mat_A,double *mat_Res,int N,int M);
extern void dgesv_( int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info );
void GetHartreePotential(double *hartree_vec,double *sij,double *kij,double *wfn,double *vec_ei,int Ne,int order,int phase)
{
	printf("HARTREE POTENTIAL MODULE\n");
	int nodes = Ne*order + 1;
	int r_nodes = nodes-2;
	int N=r_nodes;
	int NRHS=1;
	int LDA=N;
	int LDB=N;
	int ipiv[N];
	int info;
	int Norb=1;
	double sum=0.0;
	double coeffwfn,coeffrho;
	const double PI = 3.14159265358979323846;

	//LOCAL ARRAYS:
	double *rho = (double *)malloc(sizeof(double)*r_nodes);
	double *aux_kij = (double *)malloc(sizeof(double)*r_nodes*r_nodes);
	double *aux_mat = (double *)malloc(sizeof(double)*r_nodes*r_nodes);
	double *right_vec = (double *)malloc(sizeof(double)*r_nodes);
	double *aux_vec = (double *)malloc(sizeof(double)*r_nodes);
	printf("WAVE FUNCTION PHASE %d\n",phase);
	coeffwfn = 0.0;
	for(int i=0; i<r_nodes; i++)
	{
		rho[i] = 0.0;
		for(int orb=0; orb<Norb; orb++)
		{
			 rho[i] += wfn[i + orb*r_nodes]*wfn[i + orb*r_nodes];
			 //rho[i] = rho[i] + 0.5*vec_ei[orb];

		}
		rho[i] = 2.0*rho[i];
		//printf("Rho = %lf\n",rho[i]);
	}

	MatrixProduct(sij,rho,right_vec,r_nodes,r_nodes,NRHS);
	ScalarXMatrix(4.0*PI,right_vec,aux_vec,r_nodes,NRHS);
	ScalarXMatrix(2.0,kij,aux_kij,r_nodes,r_nodes);
	ColumnMayor(aux_kij,aux_mat,r_nodes,r_nodes);

	dgesv_(&N,&NRHS,aux_mat,&LDA,ipiv,aux_vec,&LDB,&info);

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
	        //printf("Hartree Potential = %lf\n",hartree_vec[i]);
	}
	free(aux_vec);
	free(right_vec);
	free(aux_kij);
	free(aux_mat);
	free(rho);

}
