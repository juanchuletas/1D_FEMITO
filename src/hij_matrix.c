#include <stdio.h>
#include <stdlib.h>
#include "data_structures.h"

void GetHmatrix(double *h,double *v,double *k,int M,int N)
{
	for(int i=0; i<M*N; i++)
	{
		h[i] = v[i] + k[i];
	}
}
