// -----------------------------------------------------------------
// Return dot product of two wilson_vectors
// Tr[adag.b]
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

complex wvec_dot(wilson_vector *a, wilson_vector *b) {
  register complex sum = cmplx(0.0, 0.0);
#ifndef FAST
  register int i;
#if DIMF > 6
  register int j;
#endif

  for(i = 0; i < 4; i++) {
    CMULJ_SUM(a->d[i].c[0], b->d[i].c[0], sum);
    CMULJ_SUM(a->d[i].c[1], b->d[i].c[1], sum);
#if DIMF > 2
    CMULJ_SUM(a->d[i].c[2], b->d[i].c[2], sum);
#if DIMF > 3
    CMULJ_SUM(a->d[i].c[3], b->d[i].c[3], sum);
#if DIMF > 4
    CMULJ_SUM(a->d[i].c[4], b->d[i].c[4], sum);
#if DIMF > 5
    CMULJ_SUM(a->d[i].c[5], b->d[i].c[5], sum);
#if DIMF > 6
    for (j = 6; j < DIMF; j++)
      CMULJ_SUM(a->d[i].c[j], b->d[i].c[j], sum);
#endif
#endif
#endif
#endif
#endif
  }

#else  // FAST version for NCOL = DIMF = 3
  register Real ar, ai, br, bi, cr, ci;

  ar = a->d[0].c[0].real;  ai = a->d[0].c[0].imag;
  br = b->d[0].c[0].real;  bi = b->d[0].c[0].imag;
  cr = ar * br + ai * bi;
  ci = ar * bi - ai * br;
  ar = a->d[0].c[1].real;  ai = a->d[0].c[1].imag;
  br = b->d[0].c[1].real;  bi = b->d[0].c[1].imag;
  cr += ar * br + ai * bi;
  ci += ar * bi - ai * br;
  ar = a->d[0].c[2].real;  ai = a->d[0].c[2].imag;
  br = b->d[0].c[2].real;  bi = b->d[0].c[2].imag;
  cr += ar * br + ai * bi;
  ci += ar * bi - ai * br;

  ar = a->d[1].c[0].real;  ai = a->d[1].c[0].imag;
  br = b->d[1].c[0].real;  bi = b->d[1].c[0].imag;
  cr += ar * br + ai * bi;
  ci += ar * bi - ai * br;
  ar = a->d[1].c[1].real;  ai = a->d[1].c[1].imag;
  br = b->d[1].c[1].real;  bi = b->d[1].c[1].imag;
  cr += ar * br + ai * bi;
  ci += ar * bi - ai * br;
  ar = a->d[1].c[2].real;  ai = a->d[1].c[2].imag;
  br = b->d[1].c[2].real;  bi = b->d[1].c[2].imag;
  cr += ar * br + ai * bi;
  ci += ar * bi - ai * br;

  ar = a->d[2].c[0].real;  ai = a->d[2].c[0].imag;
  br = b->d[2].c[0].real;  bi = b->d[2].c[0].imag;
  cr += ar * br + ai * bi;
  ci += ar * bi - ai * br;
  ar = a->d[2].c[1].real;  ai = a->d[2].c[1].imag;
  br = b->d[2].c[1].real;  bi = b->d[2].c[1].imag;
  cr += ar * br + ai * bi;
  ci += ar * bi - ai * br;
  ar = a->d[2].c[2].real;  ai = a->d[2].c[2].imag;
  br = b->d[2].c[2].real;  bi = b->d[2].c[2].imag;
  cr += ar * br + ai * bi;
  ci += ar * bi - ai * br;

  ar = a->d[3].c[0].real;  ai = a->d[3].c[0].imag;
  br = b->d[3].c[0].real;  bi = b->d[3].c[0].imag;
  cr += ar * br + ai * bi;
  ci += ar * bi - ai * br;
  ar = a->d[3].c[1].real;  ai = a->d[3].c[1].imag;
  br = b->d[3].c[1].real;  bi = b->d[3].c[1].imag;
  cr += ar * br + ai * bi;
  ci += ar * bi - ai * br;
  ar = a->d[3].c[2].real;  ai = a->d[3].c[2].imag;
  br = b->d[3].c[2].real;  bi = b->d[3].c[2].imag;
  cr += ar * br + ai * bi;
  ci += ar * bi - ai * br;

  sum.real = cr;
  sum.imag = ci;
#endif
  return sum;
}
// -----------------------------------------------------------------
