// -----------------------------------------------------------------
// Add or subtract irrep vectors
// c <-- c + b
// c <-- c - b
// c <-- a + b
// c <-- a - b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void sum_vector(vector *b, vector *c) {
  register int i;
  for (i = 0; i < DIMF; i++) {
    c->c[i].real += b->c[i].real;
    c->c[i].imag += b->c[i].imag;
  }
}

void dif_vector(vector *b, vector *c) {
  register int i;
  for (i = 0; i < DIMF; i++) {
    c->c[i].real -= b->c[i].real;
    c->c[i].imag -= b->c[i].imag;
  }
}

void add_vector(vector *a, vector *b, vector *c) {
  register int i;
  for (i = 0; i < DIMF; i++) {
    c->c[i].real = a->c[i].real + b->c[i].real;
    c->c[i].imag = a->c[i].imag + b->c[i].imag;
  }
}

void sub_vector(vector *a, vector *b, vector *c) {
  register int i;
  for (i = 0; i < DIMF; i++) {
    c->c[i].real = a->c[i].real - b->c[i].real;
    c->c[i].imag = a->c[i].imag - b->c[i].imag;
  }
}
// -----------------------------------------------------------------
