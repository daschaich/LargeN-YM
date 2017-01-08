// -----------------------------------------------------------------
// Add result of complex scalar multiplication on matrix
// c <-- a + s * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void c_scalar_mult_add_su3mat_f(su3_matrix_f *a, su3_matrix_f *b,
                                complex *s, su3_matrix_f *c) {

  register int i, j;
#ifndef NATIVEDOUBLE
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      c->e[i][j].real = a->e[i][j].real + b->e[i][j].real * s->real
                                        - b->e[i][j].imag * s->imag;
      c->e[i][j].imag = a->e[i][j].imag + b->e[i][j].imag * s->real
                                        + b->e[i][j].real * s->imag;
    }
  }
#else
  register double sr, si, br, bi, cr, ci;

  sr = (*s).real;
  si = (*s).imag;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      br = b->e[i][j].real;
      bi = b->e[i][j].imag;

      cr = sr * br - si * bi;
      ci = sr * bi + si * br;

      c->e[i][j].real = a->e[i][j].real + cr;
      c->e[i][j].imag = a->e[i][j].imag + ci;
    }
  }
#endif
}
// -----------------------------------------------------------------
