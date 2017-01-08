// -----------------------------------------------------------------
// Subtract result of complex scalar multiplication on irrep vector
// b <-- b - s * c
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void c_scalar_mult_sub_su3vec(su3_vector *b, complex *s, su3_vector *c) {
  register int i;
#ifndef NATIVEDOUBLE
  for (i = 0; i < DIMF; i++) {
    b->c[i].real -= c->c[i].real * s->real - c->c[i].imag * s->imag;
    b->c[i].imag -= c->c[i].imag * s->real + c->c[i].real * s->imag;
  }
#else
  register double sr, si, br, bi, cr, ci;

  sr = (*s).real;
  si = (*s).imag;
  for (i = 0; i < DIMF; i++) {
    br = c->c[i].real;
    bi = c->c[i].imag;

    cr = sr * br - si * bi;
    ci = sr * bi + si * br;

    b->c[i].real -= cr;
    b->c[i].imag -= ci;
  }
#endif
}
// -----------------------------------------------------------------
