// -----------------------------------------------------------------
// Complex scalar multiplication on irrep vector
// c <-- s * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void c_scalar_mult_su3vec(su3_vector *b, complex *s, su3_vector *c) {
  register int i;
#ifndef NATIVEDOUBLE
  for (i = 0; i < DIMF; i++) {
    c->c[i].real = b->c[i].real * s->real - b->c[i].imag * s->imag;
    c->c[i].imag = b->c[i].imag * s->real + b->c[i].real * s->imag;
  }
#else
  register double sr,si,br,bi,cr,ci;

  sr = (*s).real;
  si = (*s).imag;
  for (i = 0; i < DIMF; i++) {
    br=b->c[i].real;
    bi=b->c[i].imag;

    cr = sr * br - si * bi;
    ci = sr * bi + si * br;

    c->c[i].real = cr;
    c->c[i].imag = ci;
  }
#endif
}
// -----------------------------------------------------------------
