#include "data_structures.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
double Integra01( double *val){

  double v0,v1;
  double valor;

  v0 = val[0];
  v1 = val[1];

  valor = v0*v0 + v0*v1 + v1*v1;
  valor *= (2./3.);

  return valor;

}

double Integra02( double *val){

  double v0,v1,v2;
  double valor;
 
  v0 = val[0];
  v1 = val[1];
  v2 = val[2];

  valor  = 2.*v0*v0 + 2.*v0*v1 + 8.*v1*v1 ;
  valor += ( - v0*v2 + 2.*v1*v2 + 2.*v2*v2 );
  valor *= (2./15.);
  
  return valor;
}

double Integra03( double *val){

  double v0,v1,v2,v3;
  double valor;
 
  v0 = val[0];
  v1 = val[1];
  v2 = val[2];
  v3 = val[3];

  valor  = 64.*v0*v0 + 99.*v0*v1 + 324.*v1*v1 - 36.*v0*v2;
  valor += (-81.*v1*v2 + 324.*v2*v2 + 19.*v0*v3 - 36.*v1*v3);
  valor += (99.*v2*v3 + 64.*v3*v3);

  valor *= (1./420.);
  
  return valor;
}

double Integra04( double *val){

  double v0,v1,v2,v3,v4;
  double valor;
 
  v0 = val[0];
  v1 = val[1];
  v2 = val[2];
  v3 = val[3];
  v4 = val[4];

  valor  = 146.*v0*v0 + 296.*v0*v1 + 896.*v1*v1 - 174.*v0*v2;
  valor += (-384.*v1*v2 + 936.*v2*v2 + 56.*v0*v3 + 256.*v1*v2);
  valor += (-384.*v2*v3 + 896.*v3*v3 - 29.*v0*v4 + 56.*v1*v4);
  valor += (-174.*v2*v4 + 296.*v3*v4 - 146.*v4*v4);

  valor *= (2./2835.);
  
  return valor;
}

double Integra05( double *val){

  double v0,v1,v2,v3,v4,v5;
  double valor;
 
  v0 = val[0];
  v1 = val[1];
  v2 = val[2];
  v3 = val[3];
  v4 = val[4];
  v5 = val[5];

  valor  =  30512.*v0*v0 + 74325.*v0*v1 + 223250.*v1*v1 - 59550.*v0*v2;
  valor += (-147750.*v1*v2 + 249500.*v2*v2 + 34250.*v0*v3 + 114750.*v1*v3);
  valor += (-163500.*v2*v3 + 249500.*v3*v3 - 9150.*v0*v4 - 62875.*v1*v4);
  valor += (+114750.*v2*v4 - 147750.*v3*v4 + 223250.*v4*v4 + 4437.*v0*v5);
  valor += (-9150.*v1*v5 + 34250.*v2*v5 - 59550.*v3*v5 + 74325.*v4*v5);
  valor += (+30512.*v5*v5);

  valor *= (1./399168.);
  
  return valor;
}

double Integra06( double *val){

  double v0,v1,v2,v3,v4,v5,v6;
  double valor;
 
  v0 = val[0];
  v1 = val[1];
  v2 = val[2];
  v3 = val[3];
  v4 = val[4];
  v5 = val[5];
  v6 = val[6];

  valor  =  90269*v0*v0 + 252522.*v0*v1 + 776304.*v1*v1 - 254565.*v0*v2;
  valor += (-733050.*v1*v2 + 982125*v2*v2 + 204740.*v0*v3 + 682560.*v1*v3);
  valor += (-969300.*v2*v3 + 1054400*v3*v3 - 103050.*v0*v4 - 438210.*v1*v4);
  valor += (+726975.*v2*v4 - 969300.*v3*v4 + 982125.*v4*v4 + 23202.*v0*v5);
  valor += (+204768.*v1*v5 - 438219*v2*v5 + 682560.*v3*v5 - 733050.*v4*v5);
  valor += (+776304.*v5*v5 - 10237.*v0*v6 + 23202.*v1*v6 - 103050.*v2*v6);
  valor += (+204740.*v3*v6 - 254565.*v4*v6 + 252522.*v5*v6 + 90269.*v6*v6);

  valor *= (1./1501500.);
  
  return valor;
}
double IntegraElemento(int e, int p, int *g, double alpha, double *coef){

  int i;
  double phin[10];
  double valor = 0;

  for(i = 0; i < p; i++)
    phin[i] = coef[ g[p*e + i] ];


  switch(p){
    case 2: valor = Integra01(phin);
            break;

    case 3: valor = Integra02(phin);
            break;

    case 4: valor = Integra03(phin);
            break;

    case 5: valor = Integra04(phin);
            break;

    case 6: valor = Integra05(phin);
            break;

    case 7: valor = Integra06(phin);
            break;
  }

  return alpha*valor;
}
int normWF(int nele,int p,int *g, struct Element *e, double *matC){

  int i,nnt,ei,nodos,orb;
  double alpha,normC;
  double *coef;

  nnt= nele*p + 1;
  nodos = nnt - 2;
  coef = (double *)malloc(sizeof(double)*(nnt));
  for( orb = 0; orb < nodos; orb++ ){

    coef[0] = 0.;
    coef[nnt-1] = 0.;
    normC = 0;
    for( i = 1; i < nodos+1; i++)
    {
      coef[i] = matC[(i-1) + orb*nodos];
      //printf("Orb[%d],coeff[%d] = %d\n",orb,i,(i-1)*nodos + orb);
    }
    for( ei = 0;  ei < nele; ei++ ){
      alpha = 0.5*e[ei].h;
      normC += IntegraElemento(ei,p+1,g,alpha,coef);
      //printf("Orb[%d],norm = %lf\n",orb,normC);
    }

    normC = 1./sqrt(normC);
    printf("Orb[%d],norm = %lf\n",orb,normC);

    for( i = 0; i < nodos; i++)
      matC[i + orb*nodos] *= normC;

  }


  free(coef);
  return 0;
}

