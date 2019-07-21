#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "data_structures.h"
void EvaluatePotential(char kindpot[180],int Ne,int order,struct Element *e,struct Vertex *n,struct Vertex *pot)
{
	if(strcmp(kindpot,"Harmonic")==0)
	{
		switch(order)
		{
			case 1:
				for(int i=0; i<Ne; i++)
				{
					e[i].pot = (struct Vertex *)malloc(sizeof(struct Element)*(order + 1));
					for(int j=0; j<order; j++)
					{
						e[i].pot[j].x = e[i].n[j].x*e[i].n[j].x;
						e[i].pot[j+1].x = e[i].n[j+1].x*e[i].n[j+1].x;
					}
				}
				break;
			case 2:
				for(int i=0; i<Ne; i++)
				{
					e[i].pot = (struct Vertex *)malloc(sizeof(struct Element)*(order + 1));
					e[i].pot[0].x = e[i].n[0].x*e[i].n[0].x;
					e[i].pot[1].x = e[i].n[1].x*e[i].n[1].x;
					e[i].pot[2].x = e[i].n[2].x*e[i].n[2].x;

				}
		}
	}
}
