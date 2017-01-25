// -----------------------------------------------------------------
// Add result of complex scalar multiplication on irrep vector
// c <-- c + s * b
// c <-- a + s * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void c_scalar_mult_sum_su3vec(vector *b, complex *s, vector *c) {
  register int i;
  for (i = 0; i < DIMF; i++) {
    c->c[i].real += b->c[i].real * s->real - b->c[i].imag * s->imag;
    c->c[i].imag += b->c[i].imag * s->real + b->c[i].real * s->imag;
  }
}

void c_scalar_mult_add_su3vec(vector *a, vector *b, complex *s,
                              vector *c) {

  register int i;
  for (i = 0; i < DIMF; i++) {
    c->c[i].real = a->c[i].real
                 + b->c[i].real * s->real - b->c[i].imag * s->imag;
    c->c[i].imag = a->c[i].imag
                 + b->c[i].imag * s->real + b->c[i].real * s->imag;
  }
}
// -----------------------------------------------------------------
