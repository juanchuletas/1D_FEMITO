#include <stdio.h>
#include <stdlib.h>
#include "data_structures.h"

void SetElements(struct Element *e,struct Vertex *n,int Ne,double *x,int option)
{
	switch(option)
	{
		case 1:
			for(int i=0; i<Ne; i++)
			{
				e[i].n = (struct Vertex *)malloc(sizeof(struct Vertex)*2);	
		       		e[i].n[0].x = x[i];
		       		e[i].n[1].x = x[i+1];
	       			e[i].h = e[i].n[1].x - e[i].n[0].x;
        		}
	                break;
		case 2:
        		for(int i=0; i<Ne; i++)
       			{       
				e[i].n = (struct Vertex *)malloc(sizeof(struct Vertex)*3);	
				e[i].n[0].x = x[i];
			       	e[i].n[1].x = (x[i]+x[i+1])/2.0;
			       	e[i].n[2].x = x[i+1];
			      	e[i].h = e[i].n[2].x - e[i].n[0].x;
        		}
			break;



	};
}
