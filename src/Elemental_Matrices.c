#include <stdio.h>
#include <stdlib.h>
#include "Elemental_Matrices.h"
int FirstGradeElementalMatrices(int ei,int order,double coeff,int *link_mat,double *eMatS,double *eMatK,double *eMatV,double *vij)
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

	return 0;


}
int SecondGradeElementalMatrices(int ei,int order, double coeff,int *link_mat,double *eMatS, double *eMatK, double *eMatV, double *vij)
{
        int p=order+1;
        //************** Overlap elemental matrix ***************************
  	eMatS[0] = 4.0/15.0;       eMatS[1] = 2./15.0;       eMatS[2] = -1.0/15.0;
  	eMatS[3] = eMatS[1];     eMatS[4] = 16.0/15.0;      eMatS[5] = 2.0/15.0;
  	eMatS[6] = eMatS[2];     eMatS[7] = eMatS[5];     eMatS[8] = 4.0/15.0;

        //************** Kinect elemental matrix ***************************
  	eMatK[0] = 7.0/6.0;      eMatK[1] = -4.0/3.0;    eMatK[2] = 1.0/6.0;
  	eMatK[3] = eMatK[1];     eMatK[4] = 8.0/3.0;     eMatK[5] = -4.0/3.0;
  	eMatK[6] = eMatK[2];     eMatK[7] = eMatK[5];     eMatK[8] = 7.0/6.0;

        //************** Potential elemental matrix ***************************
  	eMatV[0] = (1.0/210.0)*(39.0*vij[0] + 20.0*vij[1] - 3.0*vij[2]);  
	eMatV[1] = (2.0/105.0)*( 5.0*vij[0] +  4.0*vij[1] - 2.0*vij[2]);
  	eMatV[2] = (1.0/210.0)*(-3.0*vij[0] -  8.0*vij[1] - 3.0*vij[2]);
  	eMatV[3] = eMatV[1];
  	eMatV[4] = (8.0/105.0)*(vij[0] + 12.0*vij[1] + vij[2]);
  	eMatV[5] = (-2.0/105.0)*(2.0*vij[0] - 4.0*vij[1] - 5.0*vij[2]);
  	eMatV[6] = eMatV[2];
  	eMatV[7] = eMatV[5];
  	eMatV[8] = (1./210.0)*(-3.0*vij[0] + 20.0*vij[1] + 39.0*vij[2]);

  	for(int i=0;i<p*p;i++)
	{
	       	eMatS[i] = coeff*eMatS[i];
  	        eMatK[i] = eMatK[i]/coeff;
  	        eMatV[i] = coeff*eMatV[i];
  	}
	return 0;

}
int ThirdGradeElementalMatrices(int ei,int order, double coeff,int *link_mat,double *eMatS, double *eMatK, double *eMatV, double *vij)
{
	int p=order+1;
	eMatS[0] = 16./105.;       eMatS[1] = 33./280.;    eMatS[2] = -3./70.;      eMatS[3] = 19./840.;
       	eMatS[4] = eMatS[1];       eMatS[5] = 27./35.;     eMatS[6] = -27./280.;    eMatS[7] = -3./70.;
	eMatS[8] = eMatS[2];       eMatS[9] = eMatS[6];    eMatS[10]= 27./35.;       eMatS[11]= 33./280.;
	eMatS[12]= eMatS[3];       eMatS[13]= eMatS[7];    eMatS[14]= eMatS[11];     eMatS[15]= 16./105.;


  
	eMatK[0] = 37./20.;        eMatK[1] = -189./80.;   eMatK[2] = 27./40.;    eMatK[3] = -13./80.;
	eMatK[4] = eMatK[1];       eMatK[5] = 27./5.;      eMatK[6] = -297./80.;  eMatK[7] = 27./40.;
	eMatK[8] = eMatK[2];       eMatK[9] = eMatK[6];       eMatK[10]= 27./5.;     eMatK[11]= -189./80.;
	eMatK[12]= eMatK[3];       eMatK[13]= eMatK[7];       eMatK[14]= eMatK[11];     eMatK[15]= 37./20.;

	 eMatV[0]  = (357.*vij[0] + 216.*vij[1] - 81.*vij[2] + 20.*vij[3])/3360.;
	 eMatV[1]  = (3./560.)*(12.*vij[0] + 18.*vij[1] - 9.*vij[2] + vij[3]);
	 eMatV[2]  = -3.*( 9.*vij[0] + 18.*vij[1] - 9.*vij[2] - 2.*vij[3])/1120.;
	 eMatV[3]  = (10.*vij[0] + 9.*vij[1] + 9.*vij[2] + 10.*vij[3])/1680.;
	 eMatV[4]  = eMatV[1];
	 eMatV[5]  = (27./1120.)*(4.*vij[0] + 27.*vij[1] + vij[3]);
       	 eMatV[6]  = -(27./560.)*(vij[0] + vij[3]);
       	 eMatV[7]  =  (3./1120.)*(2.*vij[0] + 9.*vij[1] - 18.*vij[2] - 9.*vij[3]);
	 eMatV[8]  = eMatV[2];
	 eMatV[9]  = eMatV[6];
       	 eMatV[10] = (27./1120.)*(vij[0] + 27.*vij[2] + 4.*vij[3]);
	 eMatV[11] = (3./560.)*(vij[0] - 9.*vij[1] + 18.*vij[2] + 12.*vij[3]);
       	 eMatV[12] = eMatV[3];
	 eMatV[13] = eMatV[7];
	 eMatV[14] = eMatV[11];
	 eMatV[15] = (20.*vij[0] - 81.*vij[1] + 216.*vij[2] + 357.*vij[3])/3360.;
	 
	 for(int i=0;i<p*p;i++)
         {
                eMatS[i] = coeff*eMatS[i];
                eMatK[i] = eMatK[i]/coeff;
                eMatV[i] = coeff*eMatV[i];
         }

	 return 0;



}

