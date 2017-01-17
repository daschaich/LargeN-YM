// -----------------------------------------------------------------
// Subtract fundamental matrices -- unlike CDIF, output is always last
// c <-- c - b
// c <-- a - b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void dif_su3_matrix_f(su3_matrix_f *b, su3_matrix_f *c) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      c->e[i][j].real -= b->e[i][j].real;
      c->e[i][j].imag -= b->e[i][j].imag;
    }
  }
}

void sub_su3_matrix_f(su3_matrix_f *a, su3_matrix_f *b, su3_matrix_f *c) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      c->e[i][j].real = a->e[i][j].real - b->e[i][j].real;
      c->e[i][j].imag = a->e[i][j].imag - b->e[i][j].imag;
    }
  }
}
// -----------------------------------------------------------------
