#include <stdio.h>
#include <stdlib.h>
#include "data_structures.h"
#include <math.h>
int GetWfnPhase(int nodos,int orb,int *phase,double *matC)
{
  int i,answer;
  double pre,val,der;
  pre = matC[orb];
  i=1;
  answer = 1;
  while( i < nodos && answer != 0 ){
    val = matC[i + orb*nodos];

    der = val - pre;

    if( fabs(der) > 1E-5 )
      answer = 0;

    pre = val;

    i++;
  }


  if( der > 0 )
    *phase = 1;
  else
    *phase = -1;


  return 0;
}
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
	fprintf(wfn_data,"Phi(r_i)\tRho(r_i)\n");
	GetWfnPhase(r_nodes,orb,&phase,wfn);
	for(int i=1; i<=r_nodes; i++)
	{
		cfwfn = wfn[(i-1) + r_nodes*orb]*phase;
		coefrho = cfwfn*cfwfn;
		cfwfn += vecE[orb];
		coefrho += vecE[orb];
		fprintf(wfn_data,"% lf\t% lf\n",cfwfn,coefrho);

		/*for(int j=0;j<order;j++)
		{
			cwfn = wfn[(k)]
			fprintf(wfn_data,"%d\t%lf\t%lf\n",k,e[i].n[j].x,wfn[k]);
			k++;

		}*/
	}
	fclose(wfn_data);
}
