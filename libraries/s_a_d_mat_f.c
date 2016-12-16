/******************  s_a_d_mat_f.c  (in su3.a) **************************
*									*
* void scalar_add_diag_su3_f(su3_matrix_f *a, Real s)	                *
* A <- A + s*I								*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void scalar_add_diag_su3_f(su3_matrix_f *a, Real s){
  register int i;

  for(i=0;i<NCOL;i++){
    a->e[i][i].real += s;
  }
}
