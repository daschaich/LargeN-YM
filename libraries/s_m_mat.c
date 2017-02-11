// -----------------------------------------------------------------
// Scalar multiplication on irrep matrix
// c <-- s * b
// c <-- c + s * b
// c <-- a + s * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void scalar_mult_mat(matrix *b, Real s, matrix *c) {
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
  c->e[0][0].real = ss * b->e[0][0].real;
  c->e[0][0].imag = ss * b->e[0][0].imag;
  c->e[0][1].real = ss * b->e[0][1].real;
  c->e[0][1].imag = ss * b->e[0][1].imag;
  c->e[0][2].real = ss * b->e[0][2].real;
  c->e[0][2].imag = ss * b->e[0][2].imag;

  c->e[1][0].real = ss * b->e[1][0].real;
  c->e[1][0].imag = ss * b->e[1][0].imag;
  c->e[1][1].real = ss * b->e[1][1].real;
  c->e[1][1].imag = ss * b->e[1][1].imag;
  c->e[1][2].real = ss * b->e[1][2].real;
  c->e[1][2].imag = ss * b->e[1][2].imag;

  c->e[2][0].real = ss * b->e[2][0].real;
  c->e[2][0].imag = ss * b->e[2][0].imag;
  c->e[2][1].real = ss * b->e[2][1].real;
  c->e[2][1].imag = ss * b->e[2][1].imag;
  c->e[2][2].real = ss * b->e[2][2].real;
  c->e[2][2].imag = ss * b->e[2][2].imag;
#endif
}

void scalar_mult_sum_mat(matrix *b, Real s, matrix *c) {
  register int i, j;
  for (i = 0; i < DIMF; i++) {
    for (j = 0; j < DIMF; j++) {
      c->e[i][j].real += s * b->e[i][j].real;
      c->e[i][j].imag += s * b->e[i][j].imag;
    }
  }
}

void scalar_mult_add_mat(matrix *a, matrix *b, Real s,
                                matrix *c) {

  register int i, j;
  for (i = 0; i < DIMF; i++) {
    for (j = 0; j < DIMF; j++) {
      c->e[i][j].real = a->e[i][j].real + s * b->e[i][j].real;
      c->e[i][j].imag = a->e[i][j].imag + s * b->e[i][j].imag;
    }
  }
}
// -----------------------------------------------------------------
