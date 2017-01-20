// -----------------------------------------------------------------
// Scalar multiplication on irrep vector
// b <- s * a
// c <- c + s * b
// c <- c - s * b
// c <- a + s * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void scalar_mult_su3_vector(su3_vector *a, Real s, su3_vector *b) {
  register int i;
  for (i = 0; i < DIMF; i++) {
    b->c[i].real = s * a->c[i].real;
    b->c[i].imag = s * a->c[i].imag;
  }
}

void scalar_mult_sum_su3_vector(su3_vector *b, Real s, su3_vector *c) {
  register int i;
  for (i = 0; i < DIMF; i++) {
    c->c[i].real += s * b->c[i].real;
    c->c[i].imag += s * b->c[i].imag;
  }
}

void scalar_mult_dif_su3_vector(su3_vector *b, Real s, su3_vector *c) {
  register int i;
  for (i = 0; i < DIMF; i++) {
    c->c[i].real -= s * b->c[i].real;
    c->c[i].imag -= s * b->c[i].imag;
  }
}

void scalar_mult_add_su3_vector(su3_vector *a, su3_vector *b, Real s,
                                su3_vector *c) {

  register int i;
  for (i = 0; i < DIMF; i++) {
    c->c[i].real = a->c[i].real + s * b->c[i].real;
    c->c[i].imag = a->c[i].imag + s * b->c[i].imag;
  }
}
// -----------------------------------------------------------------
