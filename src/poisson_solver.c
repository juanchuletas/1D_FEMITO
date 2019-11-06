#include<stdio.h>
#include<stdlib.h>
#include "data_structures.h"
#include "Elemental_Matrices.h"
extern double *ElementalOverlapMatrix(int order);
extern void MatrixProduct(double *A,double *B,double *C,int N, int M, int P);
extern void ScalarXMatrix(double coeff,double *mat_A,double *mat_Res,int N,int M);
extern void dgesv_( int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info );
void LeftHandSide(double *f,double *rho,int *link_mat,struct Element *e,int Ne,int order)
{
	int nodes=Ne*order+1;
	int p=order+1;
	int fembasis = nodes-2;
	int r_nodes = nodes-2;
	double coeff;
	double *eMatS = ElementalOverlapMatrix(order);


	for(int ei=0; ei<Ne; ei++)
	{
		coeff = 0.5*e[ei].h;
		for(int mu=0; mu<p; mu++)
		{
			f[order*ei + mu] = 0.0;
			for(int nu=0; nu<p; nu++)
			{
				int m = link_mat[p*ei+nu];
				f[order*ei + mu] += coeff*rho[m]*eMatS[mu*p+nu];
			        //printf("f[%d] = %lf\n",ei+mu,f[ei+mu]);
			}
		}
	}




}
double *PoissonSolver(double *lij,double *rho,int *link_mat,struct Element *e,int Ne,int order,int sign)
{
	int nodes = Ne*order +1;
	int fembasis_poisson = nodes-1;
	int fembasis = nodes-2;
	int NRHS = 1;
	int N=fembasis;
	int LDA=N;
	int LDB=N;
	int ipiv[N];
	int Norb=1;
	int info;
	const double PI = 3.14159265358979323846;
	const double coeff = 2.0*PI;

	
        double newrho[nodes];
        double f[nodes];
        double fnew[fembasis];
        double f_vec[fembasis];


	for(int i=0; i<fembasis; i++)
	{
		newrho[i+1] = rho[i];
	}
	newrho[0] = 0.f;
	newrho[nodes-1] = 0.f;

	LeftHandSide(f,newrho,link_mat,e,Ne,order);
	for(int i=0; i<fembasis; i++)
	{
		fnew[i] = f[i+1];
		printf("f[%d] = %lf\n",i,fnew[i]);
	}

	double *hartree_vec = (double *)malloc(sizeof(double)*fembasis);
	/*double *right_vec = (double *)malloc(sizeof(double)*fembasis_poisson);*/
	double *aux_mat = (double *)malloc(sizeof(double)*fembasis*fembasis);
	//MatrixProduct(uij,rho,right_vec,fembasis_poisson,fembasis_poisson,NRHS);
	ScalarXMatrix(2.0,lij,aux_mat,fembasis,fembasis);
	ScalarXMatrix(coeff,fnew,f_vec,fembasis,1);
	dgesv_(&N,&NRHS,aux_mat,&LDA,ipiv,fnew,&LDB,&info);

	if( info > 0 )
        {
        printf( "The diagonal element of the triangular factor of A,\n" );
        printf( "U(%i,%i) is zero, so that A is singular;\n", info, info );
        printf( "the solution could not be computed.\n" );
        exit(1);
        }

	for(int i=0; i<fembasis; i++)
	{
		hartree_vec[i] = fnew[i];
		printf("HartreePot[%d] = %lf\n",i,hartree_vec[i]);

	}
	//printf("HartreePot[%d] = %lf\n",nodes-1,hartree_vec[nodes-1]);


	//free(right_vec);
	//free(aux_mat);

	return hartree_vec;
	

}
