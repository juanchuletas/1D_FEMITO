#include <stdio.h>

void MatrixProduct(double *mat_A,double *mat_B, double *mat_C,int N,int M, int P)
{

	for(int i=0; i<N; i++)
	{
		for(int j=0; j<P; j++)
		{
			mat_C[i*P + j] = 0.0;
			for(int k=0; k<M; k++)
			{
				mat_C[i*P + j] += mat_A[i*M + k]*mat_B[j + k*P];
				
			}
		}
	}
}
void ScalarXMatrix(double coeff,double *mat_A,double *mat_Res,int N,int M)
{
	for(int i=0; i<N*M; i++)
	{
		mat_Res[i] = coeff*mat_A[i];
	}
}

