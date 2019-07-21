#include <stdio.h>
#include <stdlib.h>

#ifndef _ELEMENTAL_MATRICES_H_
#define _ELEMENTAL_MATRICES_H_
void FirstGradeElementalMatrices(int ei,int order,double coeff,int *link_mat,double *eMatS,double *eMatK,double *eMatV,double *vij)
{
        int p=order+1;
        //******* Overlap elemental matrix *****
        eMatS[0] = 2.0/3.0;  eMatS[1] = 1.0/3.0;
        eMatS[2] = eMatS[1]; eMatS[3] = eMatS[0];
        //***************************************
        //******* Kinect elemental matrix *****
        eMatK[0] = 0.5; eMatK[1] = -0.5;
        eMatK[2] = eMatK[1]; eMatK[3] = eMatK[0];
        //***************************************
        //******* Potential elemental matrix *****
        eMatV[0] = (3.0*vij[0] + vij[1])/6.0; eMatV[1] = (vij[0]+vij[1])/6.0;
        eMatV[2] = eMatV[1];  eMatV[3] = (vij[0] + vij[1]*3.0)/6.0;
        for(int i=0; i<p*p; i++)
        {
                eMatS[i] = coeff*eMatS[i];
                eMatK[i] = eMatK[i]/coeff;
                eMatV[i] = coeff*eMatV[i];
        }


}
void SecondGradeElementalMatrices(int ei,int order, double coeff,int *link_mat,double *eMatS, double *eMatK, double *eMatV, double *vij)
{
        int p=order+1;
        //************** Overlap elemental matrix ***************************
  	eMatS[0] = 4./15.;       eMatS[1] = 2./15.;       eMatS[2] = -1./15.;
  	eMatS[3] = eMatS[1];     eMatS[4] = 16./15.;      eMatS[5] = 2./15.;
  	eMatS[6] = eMatS[2];     eMatS[7] = eMatS[5];     eMatS[8] = 4./15.;

        //************** Kinect elemental matrix ***************************
  	eMatK[0] = 7./6.;     eMatK[1] = -4./3.;    eMatK[2] = 1./6.;
  	eMatK[3] = eMatK[1];     eMatK[4] = 8./3.;     eMatK[5] = -4./3.;
  	eMatK[6] = eMatK[2];     eMatK[7] = eMatK[5];     eMatK[8] = 7./6.;

        //************** Potential elemental matrix ***************************
  	eMatV[0] = (1./210.)*(39.*vij[0] + 20.*vij[1] - 3.*vij[2]);
  	eMatV[1] = (2./105.)*( 5.*vij[0] +  4.*vij[1] - 2.*vij[2]);
  	eMatV[2] = (1./210.)*(-3.*vij[0] -  8.*vij[1] - 3.*vij[2]);
  	eMatV[3] = eMatV[1];
  	eMatV[4] = (8./105.)*(vij[0] + 12.*vij[1] + vij[2]);
  	eMatV[5] = (-2./105.)*(2.*vij[0] - 4.*vij[1] - 5.*vij[2]);
  	eMatV[6] = eMatV[2];
  	eMatV[7] = eMatV[5];
  	eMatV[8] = (1./210.)*(-3.*vij[0] + 20.*vij[1] + 39.*vij[2]);

  	for(int i=0;i<p*p;i++)
	{
	       	eMatS[i] = coeff*eMatS[i];
  	        eMatK[i] = eMatK[i]/coeff;
  	        eMatV[i] = coeff*eMatV[i];
  	}

}
#endif

