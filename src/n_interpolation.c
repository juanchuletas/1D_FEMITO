#include<stdio.h>


double Li_function(int i,int N,double *x,double target)
{
	double prod=1.0;
	for(int j=0; j<N;j++)
	{
		if(j!=i)
			prod=prod*(target-x[j])/(x[i]-x[j]);

	}
	return prod;

}
double Pn_Interpolation(int order,double *x,double *y,double target)
{
	int N = order;

	double sum=0.0;
	for(int i=0;i<N;i++)
	{
		sum = sum + Li_function(i,N,x,target)*y[i];
	}

	return sum;

}
