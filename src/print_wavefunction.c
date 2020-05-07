#include <stdio.h>
#include <stdlib.h>
#include "data_structures.h"
#include "wfn_operations.h"
#include <math.h>
void PrintWaveFunction(int Ne,int order,struct Element *e,struct Vertex *n,double *wfn,double *vecE)
{
	
	FILE *wfn_data;
	double cfwfn,coefrho;
	int orb = 0;
	int nodes = Ne*order+1;
	int r_nodes = nodes-2;
	int phase=1;
	wfn_data = fopen("wfn.dat","w");
	/*if(wfn_data==NULL);
	{
		printf("ERROR IN CREATING THE WFN\n");
		exit(1);
	}*/
	printf("Wave Function in: wfn.dat\n");
	fprintf(wfn_data,"#Node\tPhi(r_i)\tRho(r_i)\n");
	GetWfnPhase(r_nodes,orb,&phase,wfn);
	printf("Wave function phase: %d\n",phase);
	for(int i=0; i<r_nodes; i++)
	{
		cfwfn = wfn[i + r_nodes*orb]*phase;
		coefrho = cfwfn*cfwfn;
		//cfwfn += 0.5*vecE[orb];
		//coefrho += 0.5*vecE[orb];
		fprintf(wfn_data,"%d\t% lf\t% lf\n",i,cfwfn,coefrho);

		/*for(int j=0;j<order;j++)
		{
			cwfn = wfn[(k)]
			fprintf(wfn_data,"%d\t%lf\t%lf\n",k,e[i].n[j].x,wfn[k]);
			k++;

		}*/
	}
	fclose(wfn_data);
}
