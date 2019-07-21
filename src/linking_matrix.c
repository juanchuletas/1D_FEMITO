#include <stdio.h>
#include <stdlib.h>

void GetLinkMatrix(int *link_mat,int Ne,int order)
{
	int p = order+1;
	for(int i=0; i<Ne; i++)
        {
                for(int j=0; j<p; j++)
                {
                        link_mat[p*i+j] = p*i+j-i;
                }
        }

}
