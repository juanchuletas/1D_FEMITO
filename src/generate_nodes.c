#include <stdio.h>
#include <stdio.h>
#include <string.h>
void GenerateNodes(double *x,int Ne, double r0, double rN, char mesh[180])
{
	double hx;
	if(strcmp(mesh,"Default\0")==0)
	{
		printf("Equally Spaced Mesh...\n");

		hx = (rN-(r0))/(double)(Ne);
		       
		for(int i=0; i<Ne+1; i++)
		{
			x[i] = r0 + hx*i;
			printf("Vertex value x[%d] = %lf\n",i,x[i]);
        
		}
	}
}
