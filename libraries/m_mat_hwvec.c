// -----------------------------------------------------------------
// Irrep matrix--half_wilson_vector multiplication
// c <-- a.b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void mult_su3_mat_hwvec(su3_matrix *a, half_wilson_vector *b,
                        half_wilson_vector *c) {
#ifndef FAST
    mult_su3_mat_vec(a, &(b->h[0]), &(c->h[0]));
    mult_su3_mat_vec(a, &(b->h[1]), &(c->h[1]));

#else   // FAST version for NCOL = DIMF = 3
  register Real a0r, a0i, a1r, a1i, a2r, a2i;
  register Real b0r, b0i, b1r, b1i, b2r, b2i;

  // h[0] multiplication
  a0r = a->e[0][0].real;    a0i = a->e[0][0].imag;
  b0r = b->h[0].c[0].real;  b0i = b->h[0].c[0].imag;
  a1r = a->e[0][1].real;    a1i = a->e[0][1].imag;
  b1r = b->h[0].c[1].real;  b1i = b->h[0].c[1].imag;
  a2r = a->e[0][2].real;    a2i = a->e[0][2].imag;
  b2r = b->h[0].c[2].real;  b2i = b->h[0].c[2].imag;
  c->h[0].c[0].real = a0r * b0r - a0i * b0i + a1r * b1r - a1i * b1i
                                            + a2r * b2r - a2i * b2i;
  c->h[0].c[0].imag = a0r * b0i + a0i * b0r + a1r * b1i + a1i * b1r
                                            + a2r * b2i + a2i * b2r;

  a0r = a->e[1][0].real;    a0i = a->e[1][0].imag;
  b0r = b->h[0].c[0].real;  b0i = b->h[0].c[0].imag;
  a1r = a->e[1][1].real;    a1i = a->e[1][1].imag;
  b1r = b->h[0].c[1].real;  b1i = b->h[0].c[1].imag;
  a2r = a->e[1][2].real;    a2i = a->e[1][2].imag;
  b2r = b->h[0].c[2].real;  b2i = b->h[0].c[2].imag;
  c->h[0].c[1].real = a0r * b0r - a0i * b0i + a1r * b1r - a1i * b1i
                                            + a2r * b2r - a2i * b2i;
  c->h[0].c[1].imag = a0r * b0i + a0i * b0r + a1r * b1i + a1i * b1r
                                            + a2r * b2i + a2i * b2r;

  a0r = a->e[2][0].real;    a0i = a->e[2][0].imag;
  b0r = b->h[0].c[0].real;  b0i = b->h[0].c[0].imag;
  a1r = a->e[2][1].real;    a1i = a->e[2][1].imag;
  b1r = b->h[0].c[1].real;  b1i = b->h[0].c[1].imag;
  a2r = a->e[2][2].real;    a2i = a->e[2][2].imag;
  b2r = b->h[0].c[2].real;  b2i = b->h[0].c[2].imag;
  c->h[0].c[2].real = a0r * b0r - a0i * b0i + a1r * b1r - a1i * b1i
                                            + a2r * b2r - a2i * b2i;
  c->h[0].c[2].imag = a0r * b0i + a0i * b0r + a1r * b1i + a1i * b1r
                                            + a2r * b2i + a2i * b2r;

  // h[1] multiplication
  a0r = a->e[0][0].real;    a0i = a->e[0][0].imag;
  b0r = b->h[1].c[0].real;  b0i = b->h[1].c[0].imag;
  a1r = a->e[0][1].real;    a1i = a->e[0][1].imag;
  b1r = b->h[1].c[1].real;  b1i = b->h[1].c[1].imag;
  a2r = a->e[0][2].real;    a2i = a->e[0][2].imag;
  b2r = b->h[1].c[2].real;  b2i = b->h[1].c[2].imag;
  c->h[1].c[0].real = a0r * b0r - a0i * b0i + a1r * b1r - a1i * b1i
                                            + a2r * b2r - a2i * b2i;
  c->h[1].c[0].imag = a0r * b0i + a0i * b0r + a1r * b1i + a1i * b1r
                                            + a2r * b2i + a2i * b2r;

  a0r = a->e[1][0].real;    a0i = a->e[1][0].imag;
  b0r = b->h[1].c[0].real;  b0i = b->h[1].c[0].imag;
  a1r = a->e[1][1].real;    a1i = a->e[1][1].imag;
  b1r = b->h[1].c[1].real;  b1i = b->h[1].c[1].imag;
  a2r = a->e[1][2].real;    a2i = a->e[1][2].imag;
  b2r = b->h[1].c[2].real;  b2i = b->h[1].c[2].imag;
  c->h[1].c[1].real = a0r * b0r - a0i * b0i + a1r * b1r - a1i * b1i
                                            + a2r * b2r - a2i * b2i;
  c->h[1].c[1].imag = a0r * b0i + a0i * b0r + a1r * b1i + a1i * b1r
                                            + a2r * b2i + a2i * b2r;

  a0r = a->e[2][0].real;    a0i = a->e[2][0].imag;
  b0r = b->h[1].c[0].real;  b0i = b->h[1].c[0].imag;
  a1r = a->e[2][1].real;    a1i = a->e[2][1].imag;
  b1r = b->h[1].c[1].real;  b1i = b->h[1].c[1].imag;
  a2r = a->e[2][2].real;    a2i = a->e[2][2].imag;
  b2r = b->h[1].c[2].real;  b2i = b->h[1].c[2].imag;
  c->h[1].c[2].real = a0r * b0r - a0i * b0i + a1r * b1r - a1i * b1i
                                            + a2r * b2r - a2i * b2i;
  c->h[1].c[2].imag = a0r * b0i + a0i * b0r + a1r * b1i + a1i * b1r
                                            + a2r * b2i + a2i * b2r;
#endif
}
// -----------------------------------------------------------------
