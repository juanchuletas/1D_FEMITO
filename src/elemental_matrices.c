#include <stdio.h>
#include <stdlib.h>


void GetElementalMatrices(int order, double *eMatS,double *eMatK, double *eMatV,double *v)
{
	switch(order)
	{
		case 1:

			//******* Overlap elemental matrix *****
			eMatS[0] = 1.0/3.0;  eMatS[1] = 1.0/6.0;
			eMatS[2] = eMatS[1]; eMatS[3] = eMatS[0];
			//***************************************
			//******* Kinect elemental matrix *****
			eMatK[0] = 1.0; eMatK[1] = -1.0;
			eMatK[2] = eMatK[1]; eMatK[3] = eMatK[0];
			//***************************************
			//******* Potential elemental matrix *****
			eMatV[0] = 0.5*v[0] + v[1]*(1.0/6.0); eMatV[1] = (1.0/6.0)*(v[0]+v[1]);
			eMatV[2] = eMatV[1];  eMatV[3] = 0.5*v[1] + (1.0/6.0)*v[0];
			break;
			
	}
	

}
