// -----------------------------------------------------------------
// Dot product of two wilson_vectors
// Return real part of dot product of two wilson_vectors
// ReTr[adag.b]
// Tr[adag.b]
// c <-- c + ReTr[adag.b]
// c <-- c + Tr[adag.b]
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

Real wvec_rdot(wilson_vector *a, wilson_vector *b) {
  register Real sum = 0.0;
#ifndef FAST
  register int i;
#if DIMF > 6
  register int j;
#endif
  register Real tr;

  for (i = 0; i < 4; i++) {
    tr = a->d[i].c[0].real * b->d[i].c[0].real; sum += tr;
    tr = a->d[i].c[0].imag * b->d[i].c[0].imag; sum += tr;
    tr = a->d[i].c[1].real * b->d[i].c[1].real; sum += tr;
    tr = a->d[i].c[1].imag * b->d[i].c[1].imag; sum += tr;
#if DIMF > 2
    tr = a->d[i].c[2].real * b->d[i].c[2].real; sum += tr;
    tr = a->d[i].c[2].imag * b->d[i].c[2].imag; sum += tr;
#if DIMF > 3
    tr = a->d[i].c[3].real * b->d[i].c[3].real; sum += tr;
    tr = a->d[i].c[3].imag * b->d[i].c[3].imag; sum += tr;
#if DIMF > 4
    tr = a->d[i].c[4].real * b->d[i].c[4].real; sum += tr;
    tr = a->d[i].c[4].imag * b->d[i].c[4].imag; sum += tr;
#if DIMF > 5
    tr = a->d[i].c[5].real * b->d[i].c[5].real; sum += tr;
    tr = a->d[i].c[5].imag * b->d[i].c[5].imag; sum += tr;
#if DIMF > 6
    for (j = 6; j < DIMF; j++) {
      tr = a->d[i].c[j].real * b->d[i].c[j].real; sum += tr;
      tr = a->d[i].c[j].imag * b->d[i].c[j].imag; sum += tr;
    }
#endif
#endif
#endif
#endif
#endif
  }
#else  // FAST version for NCOL = DIMF = 3
  register Real ar, ai, br, bi;

  ar = a->d[0].c[0].real;  ai = a->d[0].c[0].imag;
  br = b->d[0].c[0].real;  bi = b->d[0].c[0].imag;
  sum = ar * br + ai * bi;
  ar = a->d[0].c[1].real;  ai = a->d[0].c[1].imag;
  br = b->d[0].c[1].real;  bi = b->d[0].c[1].imag;
  sum += ar * br + ai * bi;
  ar = a->d[0].c[2].real;  ai = a->d[0].c[2].imag;
  br = b->d[0].c[2].real;  bi = b->d[0].c[2].imag;
  sum += ar * br + ai * bi;

  ar = a->d[1].c[0].real;  ai = a->d[1].c[0].imag;
  br = b->d[1].c[0].real;  bi = b->d[1].c[0].imag;
  sum += ar * br + ai * bi;
  ar = a->d[1].c[1].real;  ai = a->d[1].c[1].imag;
  br = b->d[1].c[1].real;  bi = b->d[1].c[1].imag;
  sum += ar * br + ai * bi;
  ar = a->d[1].c[2].real;  ai = a->d[1].c[2].imag;
  br = b->d[1].c[2].real;  bi = b->d[1].c[2].imag;
  sum += ar * br + ai * bi;

  ar = a->d[2].c[0].real;  ai = a->d[2].c[0].imag;
  br = b->d[2].c[0].real;  bi = b->d[2].c[0].imag;
  sum += ar * br + ai * bi;
  ar = a->d[2].c[1].real;  ai = a->d[2].c[1].imag;
  br = b->d[2].c[1].real;  bi = b->d[2].c[1].imag;
  sum += ar * br + ai * bi;
  ar = a->d[2].c[2].real;  ai = a->d[2].c[2].imag;
  br = b->d[2].c[2].real;  bi = b->d[2].c[2].imag;
  sum += ar * br + ai * bi;

  ar = a->d[3].c[0].real;  ai = a->d[3].c[0].imag;
  br = b->d[3].c[0].real;  bi = b->d[3].c[0].imag;
  sum += ar * br + ai * bi;
  ar = a->d[3].c[1].real;  ai = a->d[3].c[1].imag;
  br = b->d[3].c[1].real;  bi = b->d[3].c[1].imag;
  sum += ar * br + ai * bi;
  ar = a->d[3].c[2].real;  ai = a->d[3].c[2].imag;
  br = b->d[3].c[2].real;  bi = b->d[3].c[2].imag;
  sum += ar * br + ai * bi;
#endif

  return sum;
}

complex wvec_dot(wilson_vector *a, wilson_vector *b) {
  register complex sum = cmplx(0.0, 0.0);
#ifndef FAST
  register int i;
#if DIMF > 6
  register int j;
#endif

  for (i = 0; i < 4; i++) {
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

void wvec_rdot_sum(wilson_vector *a, wilson_vector *b, Real *c) {
  register int i;
#if DIMF > 6
  register int j;
#endif
  register Real tr;

  for (i = 0; i < 4; i++) {
    tr = a->d[i].c[0].real * b->d[i].c[0].real; *c += tr;
    tr = a->d[i].c[0].imag * b->d[i].c[0].imag; *c += tr;
    tr = a->d[i].c[1].real * b->d[i].c[1].real; *c += tr;
    tr = a->d[i].c[1].imag * b->d[i].c[1].imag; *c += tr;
#if DIMF > 2
    tr = a->d[i].c[2].real * b->d[i].c[2].real; *c += tr;
    tr = a->d[i].c[2].imag * b->d[i].c[2].imag; *c += tr;
#if DIMF > 3
    tr = a->d[i].c[3].real * b->d[i].c[3].real; *c += tr;
    tr = a->d[i].c[3].imag * b->d[i].c[3].imag; *c += tr;
#if DIMF > 4
    tr = a->d[i].c[4].real * b->d[i].c[4].real; *c += tr;
    tr = a->d[i].c[4].imag * b->d[i].c[4].imag; *c += tr;
#if DIMF > 5
    tr = a->d[i].c[5].real * b->d[i].c[5].real; *c += tr;
    tr = a->d[i].c[5].imag * b->d[i].c[5].imag; *c += tr;
#if DIMF > 6
    for (j = 6; j < DIMF; j++) {
      tr = a->d[i].c[j].real * b->d[i].c[j].real; *c += tr;
      tr = a->d[i].c[j].imag * b->d[i].c[j].imag; *c += tr;
    }
#endif
#endif
#endif
#endif
#endif
  }
}

void wvec_dot_sum(wilson_vector *a, wilson_vector *b, complex *c) {
  register int i;
#if DIMF > 6
  register int j;
#endif

  for (i = 0; i < 4; i++) {
    CMULJ_SUM(a->d[i].c[0], b->d[i].c[0], *c);
    CMULJ_SUM(a->d[i].c[1], b->d[i].c[1], *c);
#if DIMF > 2
    CMULJ_SUM(a->d[i].c[2], b->d[i].c[2], *c);
#if DIMF > 3
    CMULJ_SUM(a->d[i].c[3], b->d[i].c[3], *c);
#if DIMF > 4
    CMULJ_SUM(a->d[i].c[4], b->d[i].c[4], *c);
#if DIMF > 5
    CMULJ_SUM(a->d[i].c[5], b->d[i].c[5], *c);
#if DIMF > 6
    for (j = 6; j < DIMF; j++)
      CMULJ_SUM(a->d[i].c[j], b->d[i].c[j], *c);
#endif
#endif
#endif
#endif
#endif
  }
}
// -----------------------------------------------------------------
