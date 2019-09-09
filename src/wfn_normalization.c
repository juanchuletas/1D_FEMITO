#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "data_structures.h"
#include "Elemental_Matrices.h"
extern double *ElementalOverlapMatrix(int order);
double FirstGradeIntegration(int ei,int order,double *eMatS,int *link_mat,double alpha,double *cf)
{
	int p = order + 1;
	double phi[p];
	for(int i=0; i<p; i++)
		phi[i] = cf[link_mat[p*ei+i]];
	double value=0.0,v[p];
	for(int mu=0; mu<p; mu++)
	{
		for(int nu=0; nu<p; nu++)
		{
			value = value + phi[mu]*phi[nu]*(alpha*eMatS[p*mu+nu]);
		}
	}

	return value;

}
double SecondGradeIntegration(int ei,int order,int *link_mat,double alpha,double *cf)
{
        int p = order + 1;
        double phi[p];
        for(int i=0; i<p; i++)
                phi[i] = cf[link_mat[p*ei+i]];
        double value=0.0,v[p];
        double *eMatS = (double *)malloc(sizeof(double)*(p*p));
        double *eMatV = (double *)malloc(sizeof(double)*(p*p));
        double *eMatK = (double *)malloc(sizeof(double)*(p*p));
        SecondGradeElementalMatrices(ei,order,alpha,link_mat,eMatS,eMatK,eMatV,v);
        for(int mu=0; mu<p; mu++)
        {
                for(int nu=0; nu<p; nu++)
                {
                        value = value + phi[mu]*phi[nu]*eMatS[p*mu+nu];
                }
        }
        free(eMatS);
        free(eMatV);
        free(eMatK);

        return value;

}
double ThirdGradeIntegration(int ei,int order,int *link_mat,double alpha,double *cf)
{
        int p = order + 1;
        double phi[p];
        for(int i=0; i<p; i++)
                phi[i] = cf[link_mat[p*ei+i]];
        double value=0.0,v[p];
        double *eMatS = (double *)malloc(sizeof(double)*(p*p));
        double *eMatV = (double *)malloc(sizeof(double)*(p*p));
        double *eMatK = (double *)malloc(sizeof(double)*(p*p));
        ThirdGradeElementalMatrices(ei,order,alpha,link_mat,eMatS,eMatK,eMatV,v);
        for(int mu=0; mu<p; mu++)
        {
                for(int nu=0; nu<p; nu++)
                {
                        value = value + phi[mu]*phi[nu]*eMatS[p*mu+nu];
                }
        }
        free(eMatS);
        free(eMatV);
        free(eMatK);

        return value;

}
void NormWfn(double *wfn_vec,int *link_mat,struct Element *e,int order,int Ne)
{
	int p=order+1;
	int nodes = Ne*order+1; //TOTAL NODES
	int r_nodes = nodes-2; //REAL NODESa
	double Norm,ai;
	double *cf = (double *)malloc(sizeof(double)*(nodes));
	double *eMatS;
	eMatS = ElementalOverlapMatrix(order);
	cf[0]=0.0; cf[nodes-1]=0.0;
	switch(order)
	{
		case 1:
			for(int orb=0; orb<r_nodes; orb++)
        		{
                		for(int i=1; i<r_nodes+1; i++)
                		{
				      	cf[i] = wfn_vec[(i-1) + r_nodes*orb];
				}
				Norm=0.0;
				for(int ei=0; ei<Ne; ei++)
				{
					ai = 0.5*e[ei].h;
					Norm = Norm + FirstGradeIntegration(ei,order,eMatS,link_mat,ai,cf);
				}
				Norm = 1.0/sqrt(Norm);
				for(int i=0; i<r_nodes; i++)
				{
					wfn_vec[i + r_nodes*orb] *= Norm;
				}
			}
			break;
		case 2:
			for(int orb=0; orb<r_nodes; orb++)
                        {
                                for(int i=1; i<r_nodes+1; i++)
                                {
                                        cf[i] = wfn_vec[(i-1) + r_nodes*orb];
                                }
                                Norm=0.0;
                                for(int ei=0; ei<Ne; ei++)
                                {
                                        ai = 0.5*e[ei].h;
                                        Norm = Norm + SecondGradeIntegration(ei,order,link_mat,ai,cf);
                                }
                                Norm = 1.0/sqrt(Norm);
                                for(int i=0; i<r_nodes; i++)
                                {
                                        wfn_vec[i + r_nodes*orb] *= Norm;
                                }
                        }
                        break;
		case 3:
			for(int orb=0; orb<r_nodes; orb++)
                        {
                                for(int i=1; i<r_nodes+1; i++)
                                {
                                        cf[i] = wfn_vec[(i-1) + r_nodes*orb];
                                }
                                Norm=0.0;
                                for(int ei=0; ei<Ne; ei++)
                                {
                                        ai = 0.5*e[ei].h;
                                        Norm = Norm + ThirdGradeIntegration(ei,order,link_mat,ai,cf);
                                }
                                Norm = 1.0/sqrt(Norm);
                                for(int i=0; i<r_nodes; i++)
                                {
                                        wfn_vec[i + r_nodes*orb] *= Norm;
                                }
                        }
                        break;




	}
	free(eMatS);
	free(cf);

}
