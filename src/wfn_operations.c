#include <stdio.h>
#include <stdlib.h>
#include "wfn_operations.h"
#include <math.h>
int GetWfnPhase(int nodos,int orb,int *phase,double *matC)
{
  int i,answer;
  double pre,val,der;
  pre = matC[orb];
  i=1;
  answer = 1;
  while( i < nodos && answer != 0 ){
    val = matC[i + orb*nodos];

    der = val - pre;

    if( fabs(der) > 1E-5 )
      answer = 0;

    pre = val;

    i++;
  }


  if( der > 0 )
    *phase = 1;
  else
    *phase = -1;


  return 0;
}
int GetWfnRange(int nodos,int orb,int sign,double *min,double *max,double *matC,double ener){
  int i;
  double xmin,xmax,valor,valor2,coef,xmax2;
  xmin = 1.E10;
  xmax = -1.E10;
  xmax2= -1.E10;

  for(i=0;i<nodos;i++){
    coef   = matC[i*nodos + orb];
    valor  = coef + ener;
    valor2 = coef*coef + ener;
    if( valor > xmax)
      xmax = valor;
    if( valor2 > xmax2)
      xmax2 = valor2;
    if( valor < xmin)
      xmin = valor;
  }

  if( xmax < xmax2 )
    xmax = xmax2;

  valor = 0.5*fabs(xmin - xmax);

  *min = xmin - valor;
  *max = xmax + valor;


  return 0;
}

