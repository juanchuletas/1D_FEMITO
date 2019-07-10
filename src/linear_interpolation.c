#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "data_structures.h"


extern void GenerateNodes(double *x,int Ne,double r0,double rN,char mesh[180]);
extern void SetElements(struct Element *e,struct Vertex *n,int Ne, double *x,int order);

void LinearInterpolation(int Ne,double r0,double rN,char atPot[180],char mesh[180],int order)
{

	int matE = order + 1;
	double *s_mat,*x,*se_mat;

	struct Element *e;


        e = (struct Element *)malloc(sizeof(struct Element)*Ne);



	x = (double *)malloc(sizeof(double)*(Ne+1));
	s_matE = (double *)malloc(sizeof(double)*matE);

	GenerateNodes(x,Ne,r0,rN,mesh);
	SetElements(e,e->n,Ne,x,order);
	for(int i=0; i<Ne; i++)
        {
                        printf("e[%d].h = %lf\n",i,e[i].h);
        }
        for(int i=0; i<Ne; i++)
        {
                for(int j=0; j<2; j++)
                {
                        printf("e[%d].n[%d].x = %lf\n",i,j,e[i].n[j].x);

                }
        }

	/*s_mat = (double *)malloc(sizeof(double)*(Ne+1)*(Ne+1));
	se_mat = (double)malloc(sizeof(double)*())*/

	//SetZeroMat(matS,M,N);	// FILL WITH 0.0 THE MATRIX S


}
