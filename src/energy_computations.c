#include<stdio.h>
#include<stdlib.h>
#include<math.h>
void EnergyResults(double *eigenval,int N)
{
	double exact[N];
	printf("****************************************************************************\n");
	printf("State\t\tApproximate Value\tExact Value\t    Relative Error:\n");
	for(int i=0; i<N; i++)
	{
		exact[i]=0.5*(-1.0/((i+1)*(i+1)));
		eigenval[i] = 0.5*eigenval[i];
		printf("%d\t% 20.10lf\t% 20.10lf\t% 20.10E\n",i+1,eigenval[i],exact[i],fabs((exact[i]-eigenval[i])/exact[i]));
	}
	printf("****************************************************************************\n");

}
