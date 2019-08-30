#include <stdio.h>
#include <stdlib.h>

#ifndef _WFN_OP_
int GetWfnPhase(int nodos,int orb,int *phase,double *matC);
int GetWfnRange(int nodos,int orb,int sign,double *min,double *max,double *matC,double ener);
 #define _WFN_OP_
#endif
