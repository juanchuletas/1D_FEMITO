#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "data_structures.h"
extern double Pn_Interpolation(int order,double *x,double *y,double target);
extern int AtomicNumber(char target[3]);
extern double nPolExtrapolation(double *x,double *y, int N,double target);
double HarmonicPot(double this)
{
	return this*this;
}
double AtomicPot(int angular,char atom[3],double this)
{
	double r=0.12;
	double l=(double)angular;
	double z=(double)AtomicNumber(atom);
	return ((l*(l+1))/this*this - 2.0*z/this);
}
void EvaluatePotential(char kindpot[180],char atom[3],int angular,int Ne,int order,struct Element *e,struct Vertex *n,struct Vertex *pot)
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
						e[i].pot[j].x = HarmonicPot(e[i].n[j].x);
						e[i].pot[j+1].x = HarmonicPot(e[i].n[j+1].x);
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
				break;
		}
	}
	else if(strcmp(kindpot,"Atomic")==0)
	{
		int poly=7;
		double x[poly],y[poly];
		switch(order)
		{
			case 1:
				for(int i=0;i<Ne; i++)
				{
					e[i].pot = (struct Vertex *)malloc(sizeof(struct Element)*(order + 1));
					for(int j=0; j<order; j++)
					{
						e[i].pot[j].x = AtomicPot(angular,atom,e[i].n[j].x);
						e[i].pot[j+1].x = AtomicPot(angular,atom,e[i].n[j+1].x);
					}
				}
				if(e[0].n[0].x<=1.e-7)
				{
					//double x[order+1],y[order+1];
					x[0]=e[1].n[0].x; y[0]=e[1].pot[0].x;
					x[1]=e[2].n[0].x; y[1]=e[2].pot[0].x;
					x[2]=e[3].n[0].x; y[2]=e[3].pot[0].x;
					x[3]=e[4].n[0].x; y[3]=e[4].pot[0].x;
					x[4]=e[5].n[0].x; y[4]=e[5].pot[0].x;
					x[5]=e[6].n[0].x; y[5]=e[6].pot[0].x;
					x[6]=e[7].n[0].x; y[6]=e[7].pot[0].x;
					e[0].pot[0].x = Pn_Interpolation(poly,x,y,0.000000);
					printf("Potential Interpolation Value at V(0): %lf\n",e[0].pot[0].x);
				}
				break;
			case 2: 
				for(int i=0; i<Ne; i++)
				{
					e[i].pot = (struct Vertex *)malloc(sizeof(struct Element)*(order + 1));
					e[i].pot[0].x = AtomicPot(angular,atom,e[i].n[0].x);
					e[i].pot[1].x = AtomicPot(angular,atom,e[i].n[1].x);
					e[i].pot[2].x = AtomicPot(angular,atom,e[i].n[2].x);
				}
					/*double x[2],y[2];
					x[0]=e[0].n[1].x; y[0]=e[0].pot[1].x;
					x[1]=e[0].n[2].x; y[1]=e[0].pot[2].x;
					e[0].pot[0].x = Pn_Interpolation(2,x,y,0.000000);
					printf("Potential Interpolation Value at V(0): %lf\n",e[0].pot[0].x);*/
				//double x[order+1],y[order+1];
				x[0]=e[0].n[1].x; y[0]=e[0].pot[1].x;
				x[1]=e[0].n[2].x; y[1]=e[0].pot[2].x;
				x[2]=e[1].n[1].x; y[2]=e[1].pot[1].x;
				x[3]=e[1].n[2].x; y[3]=e[1].pot[2].x;
				x[4]=e[2].n[1].x; y[4]=e[2].pot[1].x;
				x[5]=e[2].n[2].x; y[5]=e[2].pot[2].x;
				x[6]=e[3].n[1].x; y[6]=e[3].pot[1].x;
				e[0].pot[0].x = Pn_Interpolation(poly,x,y,0.000000);
				printf("Interpolation Value = %lf\n",e[0].pot[0].x);
				break;
			case 3:
				for(int i=0; i<Ne; i++)
				{
					e[i].pot = (struct Vertex *)malloc(sizeof(struct Element)*(order + 1));
					e[i].pot[0].x = AtomicPot(angular,atom,e[i].n[0].x);
					e[i].pot[1].x = AtomicPot(angular,atom,e[i].n[1].x);
					e[i].pot[2].x = AtomicPot(angular,atom,e[i].n[2].x);
					e[i].pot[3].x = AtomicPot(angular,atom,e[i].n[3].x);
				}
				//double x[order+1],y[order+1];
				x[0]=e[0].n[1].x; y[0]=e[0].pot[1].x;
				x[1]=e[0].n[2].x; y[1]=e[0].pot[2].x;
				x[2]=e[0].n[3].x; y[2]=e[0].pot[3].x;
				x[3]=e[1].n[1].x; y[3]=e[1].pot[1].x;
				x[4]=e[1].n[2].x; y[4]=e[1].pot[2].x;
				x[5]=e[1].n[3].x; y[5]=e[1].pot[3].x;
				x[6]=e[2].n[1].x; y[6]=e[2].pot[1].x;
				//e[0].pot[0].x = Pn_Interpolation(poly,x,y,0.000000);
				e[0].pot[0].x = nPolExtrapolation(x,y,poly,0.000000);
				printf("Potential Interpolation Value V(0): %lf\n",e[0].pot[0].x);
				//printf("Interpolation Value = %lf\n",e[0].pot[0].x);
				break;
		}
	}
}
