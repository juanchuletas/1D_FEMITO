#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern void LagrangeInterpolation(int Ne,double r0,double rN,char atPot[180],char mesh[180],char atom[3],int order,int angular);
void load_data(int Ne,double r0,double rN,char atPot[180],int order,char mesh[180],char atom[3],int angular)
{
	LagrangeInterpolation(Ne,r0,rN,atPot,mesh,atom,order,angular);
}
int main(int argc, char **argv)
{
	printf("****************************************************************************\n");
        printf("                      ********** FEMITO **********                          \n");
        printf("****************************************************************************\n");
        printf("********            Hartree-Fock Approximation In The               ********\n");
        printf("********                 Finite Element Method                      ********\n");
        printf("********                        For Atoms                           ********\n");
        printf("********          With Lagrange Interpolation Functions             ********\n");
	printf("****************************************************************************\n");
	printf("**********************     Author: Juan Garcia     *************************\n");
        printf("****************************************************************************\n");
	FILE *input;
	char fname[1024],str1[180],str2[180],str3[180];
	char atPot[180],mesh[180],atom[3];
	char *DatPot[]={"Harmonic","Atomic"};
	char *DatMesh[]={"Uniform","Froese-Fischer"};
	int Ne,order;
	double r0,rN;
	int angular=-1;
	//printf("Enter the input file with data:\n");
	//scanf("%s",fname);


	input = fopen(argv[1],"r");
	if(input==NULL)
	{
		perror("FILE NOT FOUND");
		exit(1);
	}
	//Potential: char <Harmonic,Atomic>
	fscanf(input,"%s %s",str1,atPot);
	if(strcmp(atPot,DatPot[0])!=0 && strcmp(atPot,DatPot[1])!=0)
	{
		printf("ERROR IN POTENTIAL INPUT LINE\n");
		printf("AVAILABLE POTENTIALS:\n");
		printf("%s\n",DatPot[0]);
		printf("%s\n",DatPot[1]);
		exit(1);
	}
	printf("Potential: %s\n",atPot);
	if(strcmp(atPot,"Atomic")==0)
	{
		//Atom: char <H,He>
		fscanf(input,"%s %s",str1,atom);
		printf("Atom: %s\n",atom);
		//Angular Momentum (l): char <H,He>
		fscanf(input,"%s %s %s %d",str1,str2,str3,&angular);
		printf("Angular Momentum: %d\n",angular);

	}
	//Number of Elements: int
	fscanf(input,"%s %s %s %d",str1,str2,str3,&Ne);
	printf("Number Of Elements: %d\n",Ne);
	//Initial Point: int
	fscanf(input,"%s %s %lf",str1,str2,&r0);
	//Final Point: int
	fscanf(input,"%s %s %lf",str1,str2,&rN);
	printf("Initial And Final Points: [%lf,%lf]\n",r0,rN);
	//Lagrange Interpolation Order: int <1,2>
	fscanf(input,"%s %s %s %d",str1,str2,str3,&order);
	printf("Lagrange Interpolation Order: %d\n",order);
	//Kind Of Mesh: <Default,Frose-Fischer> 
	fscanf(input,"%s %s %s %s",str1,str2,str3,mesh);
	if(strcmp(mesh,DatMesh[0])!=0 && strcmp(mesh,DatMesh[1])!=0)
	{
		printf("ERROR IN MESH INPUT LINE\n");
		printf("AVAILABLE POTENTIALS:\n");
		printf("%s\n",DatMesh[0]);
		printf("%s\n",DatMesh[1]);
		exit(1);
	}
	printf("Kind Of Mesh: %s\n",mesh);

	fclose(input);

	load_data(Ne,r0,rN,atPot,order,mesh,atom,angular);
	return 0;

}
