#include <stdio.h>
#include <stdlib.h>

#ifndef _ELEMENTAL_MAT_
 #define _ELEMENTAL_MAT_
int FirstGradeElementalMatrices(int ei,int order,double coeff,int *link_mat,double *eMatS,double *eMatK,double *eMatV,double *vij);
int SecondGradeElementalMatrices(int ei,int order, double coeff,int *link_mat,double *eMatS, double *eMatK, double *eMatV, double *vij);
int ThirdGradeElementalMatrices(int ei,int order, double coeff,int *link_mat,double *eMatS, double *eMatK, double *eMatV, double *vij);
#endif
