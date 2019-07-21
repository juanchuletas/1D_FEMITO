#include <stdio.h>
#include <stdlib.h>

extern void LagrangeInterpolation(int Ne,double r0,double rN,char atPot[180],char mesh[180],int order);
void load_data(int Ne,double r0,double rN,char atPot[180],int order,char mesh[180])
{
	LagrangeInterpolation(Ne,r0,rN,atPot,mesh,order);
}
int main(int argc, char **argv)
{
	printf("*********************************************************************\n");
        printf("                   ********** FEMITO **********                      \n");
        printf("*********************************************************************\n");
        printf("*****            Hartree-Fock Approximation In The              *****\n");
        printf("*****                 Finite Element Method                     *****\n");
        printf("*****                        For Atoms                          *****\n");
        printf("*****          With Lagrange Interpolation Functions            *****\n");
        printf("*********************************************************************\n");
	FILE *input;
	char fname[1024],str1[180],str2[180],str3[180];
	char atPot[180],mesh[180];
	int Ne,order;
	double r0,rN;
	//printf("Enter the input file with data:\n");
	//scanf("%s",fname);


	input = fopen(argv[1],"r");
	if(input==NULL)
	{
		perror("FILE NOT FOUND");
		exit(1);
	}
	fscanf(input,"%s %s %d",str1,str2,&Ne);
	fscanf(input,"%s %s %lf",str1,str2,&r0);
	fscanf(input,"%s %s %lf",str1,str2,&rN);
	fscanf(input,"%s %s",str1,atPot);
	fscanf(input,"%s %s %s %d",str1,str2,str3,&order);
	fscanf(input,"%s %s %s %s",str1,str2,str3,mesh);

	fclose(input);

	load_data(Ne,r0,rN,atPot,order,mesh);
	return 0;

}
