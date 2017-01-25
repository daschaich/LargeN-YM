// -----------------------------------------------------------------
// Traces of two fundamental matrices
// Return ReTr[adag.b]
// c <-- c + ReTr[adag.b]
// Return Tr[adag.b]
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

Real realtrace_f(matrix_f *a, matrix_f *b) {
  register int i, j;
  register Real sum = 0.0;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      sum += a->e[i][j].real * b->e[i][j].real;
      sum += a->e[i][j].imag * b->e[i][j].imag;
    }
  }
  return sum;
}

void realtrace_sum_f(matrix_f *a, matrix_f *b, Real *c) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      *c += a->e[i][j].real * b->e[i][j].real;
      *c += a->e[i][j].imag * b->e[i][j].imag;
    }
  }
}

complex complextrace_f(matrix_f *a, matrix_f *b) {
  register int i,j;
  complex sum = cmplx(0.0, 0.0);
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      sum.real += a->e[i][j].real * b->e[i][j].real;
      sum.real += a->e[i][j].imag * b->e[i][j].imag;
      sum.imag += a->e[i][j].real * b->e[i][j].imag;
      sum.imag -= a->e[i][j].imag * b->e[i][j].real;
    }
  }
  return sum;
}
// -----------------------------------------------------------------
