#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void EvaluatePotential(char pot[180],int nodes,double *x,double *v)
{
	if(strcmp(pot,"Harmonic")==0)
	{
		for(int i=0; i<nodes; i++)
		{
			v[i] = x[i]*x[i];
			printf("V(%d) = %lf\n",i,v[i]);
		}
	}
}
