/******************  dumpmat_f.c  (in su3.a) ******************************
*									*
*  void dumpmat_f( su3_matrix_f *mat )					*
*  print out the anti-hermitian part of an NCOLxNCOL complex matrix					*
*/
#include "../include/config.h"
#include <stdio.h>
#include "../include/complex.h"
#include "../include/su3.h"

void dump_ahmat_f( su3_matrix_f *m ){
  int i,j;
  for(i=0;i<NCOL;i++){
    for(j=0;j<NCOL;j++){
      printf("(%.8e,%.8e)   ",
          (m->e[i][j].real - m->e[j][i].real)/2.,
          (m->e[i][j].imag + m->e[j][i].imag)/2. );
    }
    printf("\n");
  }
  printf("\n");
}
