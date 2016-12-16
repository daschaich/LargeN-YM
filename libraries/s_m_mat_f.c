/******************  s_m_mat_f.c  (in su3.a) ******************************
*									*
* void scalar_mult_su3_matrix_f( su3_matrix_f *a, Real s, su3_matrix_f *b)	*
* B <- s*A								*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

/* b <- s*a, matrices */
void scalar_mult_su3_matrix_f( su3_matrix_f *a, Real s, su3_matrix_f *b ){

register int i,j;
    for(i=0;i<NCOL;i++)for(j=0;j<NCOL;j++){
	b->e[i][j].real = s*a->e[i][j].real;
	b->e[i][j].imag = s*a->e[i][j].imag;
    }
}
