#include <stdio.h>
#include <stdlib.h>


void ElementalMatrix()
{

	eMatS[0] = 1.0/3.0;  eMatS[1] = 1.0/6.0;
	eMatS[2] = eMats[1]; eMAtS[3] = eMatS[0];

	eMatK[0] = 1.0; eMat[1] = -1.0;
	eMatK[2] = eMat[1]; eMatK[3] = eMatK[0];
	
	for(int i=0; i<Ne; i++)
	{
		for(int j=0; j<eMatSize; j++)
		{
			eMatS[i] = eMatS[i]*e[j].h;
			eMatK[i] = eMatK[i]/e[j].h;
		}
		eMatV[0] = 
	}

}
