// -----------------------------------------------------------------
// Scalar multiplication on irrep matrix
// c <-- s * b
// c <-- c + s * b
// c <-- a + s * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void scalar_mult_su3_matrix(su3_matrix *b, Real s, su3_matrix *c) {
#ifndef FAST
  register int i, j;
  for (i = 0; i < DIMF; i++) {
    for (j = 0; j < DIMF; j++) {
      c->e[i][j].real = s * b->e[i][j].real;
      c->e[i][j].imag = s * b->e[i][j].imag;
    }
  }

#else  // FAST version for NCOL = DIMF = 3
  register Real ss = s;
  b->e[0][0].real = ss * a->e[0][0].real;
  b->e[0][0].imag = ss * a->e[0][0].imag;
  b->e[0][1].real = ss * a->e[0][1].real;
  b->e[0][1].imag = ss * a->e[0][1].imag;
  b->e[0][2].real = ss * a->e[0][2].real;
  b->e[0][2].imag = ss * a->e[0][2].imag;

  b->e[1][0].real = ss * a->e[1][0].real;
  b->e[1][0].imag = ss * a->e[1][0].imag;
  b->e[1][1].real = ss * a->e[1][1].real;
  b->e[1][1].imag = ss * a->e[1][1].imag;
  b->e[1][2].real = ss * a->e[1][2].real;
  b->e[1][2].imag = ss * a->e[1][2].imag;

  b->e[2][0].real = ss * a->e[2][0].real;
  b->e[2][0].imag = ss * a->e[2][0].imag;
  b->e[2][1].real = ss * a->e[2][1].real;
  b->e[2][1].imag = ss * a->e[2][1].imag;
  b->e[2][2].real = ss * a->e[2][2].real;
  b->e[2][2].imag = ss * a->e[2][2].imag;
#endif
}

void scalar_mult_sum_su3_matrix(su3_matrix *b, Real s, su3_matrix *c) {
  register int i, j;
  for (i = 0; i < DIMF; i++) {
    for (j = 0; j < DIMF; j++) {
      c->e[i][j].real += s * b->e[i][j].real;
      c->e[i][j].imag += s * b->e[i][j].imag;
    }
  }
}

void scalar_mult_add_su3_matrix(su3_matrix *a, su3_matrix *b, Real s,
                                su3_matrix *c) {

  register int i, j;
  for (i = 0; i < DIMF; i++) {
    for (j = 0; j < DIMF; j++) {
      c->e[i][j].real = a->e[i][j].real + s * b->e[i][j].real;
      c->e[i][j].imag = a->e[i][j].imag + s * b->e[i][j].imag;
    }
  }
}
// -----------------------------------------------------------------
