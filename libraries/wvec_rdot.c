// -----------------------------------------------------------------
// Return real part of dot product of two wilson_vectors
// ReTr[adag.b]
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
// -----------------------------------------------------------------
