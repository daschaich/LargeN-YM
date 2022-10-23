// -----------------------------------------------------------------
// Scalar multiplication on matrix
// c <-- s * b
// c <-- c + s * b
// c <-- c - s * b
// c <-- a + s * b
// c <-- a - s * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void scalar_mult_mat(matrix *b, Real s, matrix *c) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      c->e[i][j].real = s * b->e[i][j].real;
      c->e[i][j].imag = s * b->e[i][j].imag;
    }
  }
}

void scalar_mult_sum_mat(matrix *b, Real s, matrix *c) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      c->e[i][j].real += s * b->e[i][j].real;
      c->e[i][j].imag += s * b->e[i][j].imag;
    }
  }
}

void scalar_mult_dif_mat(matrix *b, Real s, matrix *c) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      c->e[i][j].real -= s * b->e[i][j].real;
      c->e[i][j].imag -= s * b->e[i][j].imag;
    }
  }
}

void scalar_mult_add_mat(matrix *a, matrix *b, Real s, matrix *c) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      c->e[i][j].real = a->e[i][j].real + s * b->e[i][j].real;
      c->e[i][j].imag = a->e[i][j].imag + s * b->e[i][j].imag;
    }
  }
}

void scalar_mult_sub_mat(matrix *a, matrix *b, Real s, matrix *c) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      c->e[i][j].real = a->e[i][j].real - s * b->e[i][j].real;
      c->e[i][j].imag = a->e[i][j].imag - s * b->e[i][j].imag;
    }
  }
}
// -----------------------------------------------------------------
