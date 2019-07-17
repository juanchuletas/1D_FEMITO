#include <stdio.h>
#include <stdlib.h>
#include <math.h>
extern void dsyev_(char*, char*, int*, double*, int*, double*, double*, int*, int*);
void diag (int n, double *h, double *s, double *e, double *v)
{
  int lwork, i, j, k, nn, ij, info;
  static double small = 1.0e-10;
  static char *V = "V";
  static char *U = "U";
  double *work, *aux, *uno;

  lwork=3*n;
  work = (double *) malloc( lwork * sizeof (double));

  aux = (double *) malloc( n * n * sizeof (double));
  for ( ij = 0; ij < n*n; ij++ )
  {
    aux[ij] = s[ij];
  }

  dsyev_( V, U, &n, aux, &n, e, work, &lwork, &info );

  if ( info !=0 )
  {
    fprintf (stderr, "Overlap matrix diagonalization failed\n");
    exit (1);
  }
// print_matrix( "S matrix ", nbasis, nbasis, aux, nbasis );
 //printf("\n");

 nn = 0;
  for ( i = 0; i < n; i++ )
  {
    if ( e[i] > small)
    {
      for ( j = 0; j < n; j++ )
      {
        aux[j+nn*n] = aux[j+i*n] / sqrt(e[i]);
      }
      nn++;
    }
  }
  if ( nn < n )
  {
     fprintf (stdout, " # of linearly independent vectors = %d\n", nn);
  }

 //print_matrix( "X matrix ", nbasis, nbasis, aux, nbasis );
 //printf("\n");

  for ( i = 0; i < n; i++)
  {
    for ( j = 0; j < nn; j++)
    {
      v[i+j*n] = 0.0;
      for  ( k = 0; k < n; k++)
      {
        v[i+j*n] = v[i+j*n] + h[i+k*n] * aux[k+j*n];
      }
    }
  }
  uno = (double *) malloc( nn * nn * sizeof (double));
  for ( i = 0; i < nn; i++ )
  {
    for ( j = 0; j < nn; j++ )
    {
      uno[i+j*nn] = 0.0;
      for  ( k = 0; k < n; k++)
      {
        uno[i+j*nn] = uno[i+j*nn] + aux[k+i*n] * v[k+j*n];
      }
    }
  }
  //print_matrix( "F\' matrix ", nbasis, nbasis, uno, nbasis );
  //printf("\n");
   /* ---------------- Diagonalization of transformed F or H matrix -------------*/
  dsyev_("V", "U", &nn, uno, &nn, e, work, &lwork, &info );

  if ( info !=0 )
  {
    fprintf (stderr, "Fock matrix diagonalization failed\n");
    exit (1);
  }


  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < nn; j++ )
    {
      v[i+j*n] = 0.0;
      for  ( k = 0; k < n; k++ )
      {
        v[i+j*n] = v[i+j*n] + aux[i+k*n] * uno[k+j*nn];
      }
    }
  }
	//print_matrix( "eigenvalues ", 1, nbasis, e, 1 );
  	//printf("\n");
  free(uno); free(aux); free(work);
  return;
}
