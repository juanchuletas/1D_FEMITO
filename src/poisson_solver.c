double *PoissonSolver(double *globalSmat,double *globalKmat,double *wfn,int Ne,int order)
{
	int nodes = Ne*order +1;
	int fembasis_poisson = nodes-1;

	double *hartree_vec = (double *)malloc(sizeof(double)*nodes);

	double *k_mat = (double *)malloc(sizeof(double)*(fembasis_poisson*fembasis_poisson));
	double *s_mat = (double *)malloc(sizeof(double)*(fembasis_poisson*fembasis_poisson));
	

	m = nodes;
	l=0;
  for(i=0;i<m-1;i++)
    for(j=0;j<m-1;j++){
      matSm[l] = matS[i*m+j];
      matVm[l] = matV[i*m+j];
      matKm[l] = matK[i*m+j];
      l++;
    }


}
