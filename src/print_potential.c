#include <stdio.h>
#include <stdlib.h>
#include "data_structures.h"

void PrintPotentialData(int order,int Ne, struct Element *e,struct Vertex *n,struct Vertex *pot)
{
	FILE *out;

	out=fopen("potential.dat","w");
	if(out==NULL)
	{
		printf("ERROR\n");
		exit(1);
	}
	fprintf(out,"Node\tr_i\t \tV(r_i)\n");
	int k=0;
	for(int i=0; i<Ne; i++)
	{
		for(int j=0; j<order; j++)
		{
			fprintf(out,"%d\t%lf\t%lf\n",k,e[i].n[j].x,e[i].pot[j].x);
			k++;
		}
	}
	fclose(out);
}
