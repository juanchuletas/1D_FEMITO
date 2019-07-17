#include <stdio.h>
#include <stdlib.h>


void GetElementalMatrices(int order, double *eMatS,double *eMatK, double *eMatV,double *v,double *link_mat)
{
	double vi[order+1];
	for(int i=0; i<p; i++)
		vi[i] = v[link_mat[p*e + i ] ];


	switch(order)
	{
		case 1:

			//******* Overlap elemental matrix *****
			eMatS[0] = 2.0/3.0;  eMatS[1] = 1.0/3.0;
			eMatS[2] = eMatS[1]; eMatS[3] = eMatS[0];
			//***************************************
			//******* Kinect elemental matrix *****
			eMatK[0] = 0.5; eMatK[1] = -0.5;
			eMatK[2] = eMatK[1]; eMatK[3] = eMatK[0];
			//***************************************
			//******* Potential elemental matrix *****
			eMatV[0] = (3.0*v[0] + v[1])/6.0; eMatV[1] = (v[0]+v[1])/6.0;
			eMatV[2] = eMatV[1];  eMatV[3] = (v[0] + v[1]*3.0)/6.0;
			break;
			
	}
	

}
