#include <stdio.h>
#include <stdlib.h>

void print_matrix( char* desc, int m, int n, double* a)
{
        int i, j;
        printf( "*****   %s   *****\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( "%lf\t", a[i*n+j] );
                printf( "\n" );
        }
}
