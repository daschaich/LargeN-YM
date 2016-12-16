/******************  dumpmat_f.c  (in su3.a) ******************************
*									*
*  void dumpmat_f( su3_matrix_f *mat )					*
*  print out a 3x3 complex matrix					*
*/
#include "../include/config.h"
#include <stdio.h>
#include "../include/complex.h"
#include "../include/su3.h"

void dumpmat_f( su3_matrix_f *m ){
int i,j;
    for(i=0;i<NCOL;i++){
	for(j=0;j<NCOL;j++)printf("(%.8e,%.8e)\t",
	    m->e[i][j].real,m->e[i][j].imag);
	printf("\n");
    }
    printf("\n");
}
