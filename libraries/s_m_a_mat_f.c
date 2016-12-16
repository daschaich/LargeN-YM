/****************  s_m_a_mat_f.c  (in su3.a) ******************************
*									*
* void scalar_mult_add_su3_matrix_f( su3_matrix_f *a, su3_matrix_f *b,	*
*	Real s, su3_matrix_f *c)						*
* C <- A + s*B								*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

/* c <- a + s*b, matrices */
void scalar_mult_add_su3_matrix_f(su3_matrix_f *a,su3_matrix_f *b,Real s,
	su3_matrix_f *c){

#ifndef NATIVEDOUBLE
register int i,j;
    for(i=0;i<NCOL;i++)for(j=0;j<NCOL;j++){
	c->e[i][j].real = a->e[i][j].real + s*b->e[i][j].real;
	c->e[i][j].imag = a->e[i][j].imag + s*b->e[i][j].imag;
    }

#else /* RS6000 version */

  register double ss;

  ss = s;

  c->e[0][0].real = a->e[0][0].real + ss*b->e[0][0].real;
  c->e[0][0].imag = a->e[0][0].imag + ss*b->e[0][0].imag;
  c->e[0][1].real = a->e[0][1].real + ss*b->e[0][1].real;
  c->e[0][1].imag = a->e[0][1].imag + ss*b->e[0][1].imag;
  c->e[0][2].real = a->e[0][2].real + ss*b->e[0][2].real;
  c->e[0][2].imag = a->e[0][2].imag + ss*b->e[0][2].imag;

  c->e[1][0].real = a->e[1][0].real + ss*b->e[1][0].real;
  c->e[1][0].imag = a->e[1][0].imag + ss*b->e[1][0].imag;
  c->e[1][1].real = a->e[1][1].real + ss*b->e[1][1].real;
  c->e[1][1].imag = a->e[1][1].imag + ss*b->e[1][1].imag;
  c->e[1][2].real = a->e[1][2].real + ss*b->e[1][2].real;
  c->e[1][2].imag = a->e[1][2].imag + ss*b->e[1][2].imag;

  c->e[2][0].real = a->e[2][0].real + ss*b->e[2][0].real;
  c->e[2][0].imag = a->e[2][0].imag + ss*b->e[2][0].imag;
  c->e[2][1].real = a->e[2][1].real + ss*b->e[2][1].real;
  c->e[2][1].imag = a->e[2][1].imag + ss*b->e[2][1].imag;
  c->e[2][2].real = a->e[2][2].real + ss*b->e[2][2].real;
  c->e[2][2].imag = a->e[2][2].imag + ss*b->e[2][2].imag;

#endif
}
