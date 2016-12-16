/******************  dumpmat.c  (in su3.a) ******************************
*									*
*  void dumpmat( su3_matrix *mat )					*
*  print out a 3x3 complex matrix					*
*/
#include "../include/config.h"
#include <stdio.h>
#include "../include/complex.h"
#include "../include/su3.h"

void dumpmat( su3_matrix *m ){
int i,j;
    for(i=0;i<DIMF;i++){
	for(j=0;j<DIMF;j++)printf("(%e,%e)\t",
	    m->e[i][j].real,m->e[i][j].imag);
	printf("\n");
    }
    printf("\n");
}
