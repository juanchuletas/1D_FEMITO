#include<stdio.h>
#include<stdlib.h>
extern void MatrixProduct(double *A,double *B,double *C,int N, int M, int P);
extern void ScalarXMatrix(double coeff,double *mat_A,double *mat_Res,int N,int M);
extern void dgesv_( int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info );
double *PoissonSolver(double *uij,double *lij,double *wfn,int Ne,int order,int sign)
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
	double *right_vec = (double *)malloc(sizeof(double)*fembasis);
	double *dummy = (double *)malloc(sizeof(double)*fembasis_poisson);
	double *aux_vec = (double *)malloc(sizeof(double)*fembasis_poisson);
	double *rho = (double *)malloc(sizeof(double)*fembasis);
	dummy[0] = 0.f;
	//rho[0] = 0.f;
	double coeffwfn;
	for(int i=0; i<fembasis; i++)
	{
		rho[i] = 0.0;
		for(int orb=0; orb<Norb; orb++)
		{
			coeffwfn = wfn[i + orb*fembasis]*sign;
			rho[i] += 2.0*(coeffwfn*coeffwfn);
		}
		printf("Rho[%d] = %lf\n",i,rho[i]);

	}

	MatrixProduct(uij,rho,right_vec,fembasis,fembasis,NRHS);
	//printf("Dummy[%d] = %lf\n",0,dummy[0]);
	for(int i=0; i<fembasis; i++)
	{
		dummy[i+1] = right_vec[i];
		//printf("Dummy[%d] = %lf\n",i+1,dummy[i+1]);
	}
	ScalarXMatrix(2.0,dummy,aux_vec,fembasis_poisson,NRHS);
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



	return hartree_vec;
	

	free(right_vec);
	free(rho);
	free(aux_vec);

}
