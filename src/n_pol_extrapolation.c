#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double nPolExtrapolation(double *x,double *y, int N,double target)
{
	double *c_vec,*d_vec,dy,ho,hp,den,w;
	double result,deltaTot,delta;
	dy = 0.00000000000001;
	int ns=0;
	c_vec = (double *)malloc(sizeof(double)*N);
	d_vec = (double *)malloc(sizeof(double)*N);
	delta = fabs(target-x[0]);
	for(int i=0; i<N; i++)
	{
		if((deltaTot=fabs(target-x[i])) < delta)
		{
			ns = i;
			delta = deltaTot;
		}
		c_vec[i] = y[i];
		d_vec[i] = y[i];
	}
	result = y[ns];
	ns = ns - 1;
	for(int m=0; m<N-1; m++)
	{
		for(int i=0; i<N-m-1; i++)
		{
			ho = x[i]-target;
			hp = x[i+m+1]-target;
			w = c_vec[i+1]-d_vec[i];
			if((den=ho-hp)==0.0){printf("Error in Extrapolation\n");exit(1);}
			den=w/den;
			d_vec[i] = hp*den;
			c_vec[i] = ho*den;
		}
		if( 2*(ns+1) < N-m-1 ) 
		{
			dy = c_vec[ns+1];
		} 
		else 
		{
			dy = d_vec[ns];
			ns = ns-1;
		}
		result = result + (dy);
	}
	return result;



}
