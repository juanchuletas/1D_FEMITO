#include<stdio.h>
#include<stdlib.h>
extern void MatrixProduct(double *A,double *B,double *C,int N, int M, int P);
extern void ScalarXMatrix(double coeff,double *mat_A,double *mat_Res,int N,int M);
extern void dgesv_( int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info );
double *PoissonSolver(double *uij,double *lij,double *rho,int Ne,int order,int sign)
{
	int nodes = Ne*order +1;
	int fembasis_poisson = nodes-1;
	int fembasis = nodes-2;
	int NRHS = 1;
	int N=fembasis_poisson;
	int LDA=N;
	int LDB=N;
	int ipiv[N];
	int Norb=1;
	int info;
	const double PI = 3.14159265358979323846;

	

	double *hartree_vec = (double *)malloc(sizeof(double)*nodes);
	double *right_vec = (double *)malloc(sizeof(double)*fembasis_poisson);
	double *aux_vec = (double *)malloc(sizeof(double)*fembasis_poisson);
	MatrixProduct(uij,rho,right_vec,fembasis_poisson,fembasis_poisson,NRHS);
	ScalarXMatrix(2.0,right_vec,aux_vec,fembasis_poisson,NRHS);
	dgesv_(&N,&NRHS,lij,&LDA,ipiv,aux_vec,&LDB,&info);

	if( info > 0 )
        {
        printf( "The diagonal element of the triangular factor of A,\n" );
        printf( "U(%i,%i) is zero, so that A is singular;\n", info, info );
        printf( "the solution could not be computed.\n" );
        exit(1);
        }

	for(int i=0; i<fembasis_poisson; i++)
	{
		hartree_vec[i] = aux_vec[i];
		printf("HartreePot[%d] = %lf\n",i,hartree_vec[i]);

	}
	hartree_vec[nodes-1] = 0.f;
	printf("HartreePot[%d] = %lf\n",nodes-1,hartree_vec[nodes-1]);


	free(right_vec);
	free(aux_vec);

	return hartree_vec;
	

}
