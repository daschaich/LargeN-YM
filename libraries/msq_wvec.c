// -----------------------------------------------------------------
// Squared magnitude of wilson_vector
// Return ReTr[adag.a]
// c <-- c + ReTr[adag.a]
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

Real magsq_wvec(wilson_vector *a) {
  register Real sum;
#ifndef FAST
  register int i;
  sum = magsq_su3vec(&(a->d[0]));
  for (i = 1; i < 4; i++)
    sum += magsq_su3vec(&(a->d[i]));

#else  // FAST version for NCOL = DIMF = 3
  register Real ar, ai;

  ar = a->d[0].c[0].real; ai = a->d[0].c[0].imag;
  sum = ar * ar + ai * ai;
  ar = a->d[0].c[1].real; ai = a->d[0].c[1].imag;
  sum += ar * ar + ai * ai;
  ar = a->d[0].c[2].real; ai = a->d[0].c[2].imag;
  sum += ar * ar + ai * ai;

  ar = a->d[1].c[0].real; ai = a->d[1].c[0].imag;
  sum += ar * ar + ai * ai;
  ar = a->d[1].c[1].real; ai = a->d[1].c[1].imag;
  sum += ar * ar + ai * ai;
  ar = a->d[1].c[2].real; ai = a->d[1].c[2].imag;
  sum += ar * ar + ai * ai;

  ar = a->d[2].c[0].real; ai = a->d[2].c[0].imag;
  sum += ar * ar + ai * ai;
  ar = a->d[2].c[1].real; ai = a->d[2].c[1].imag;
  sum += ar * ar + ai * ai;
  ar = a->d[2].c[2].real; ai = a->d[2].c[2].imag;
  sum += ar * ar + ai * ai;

  ar = a->d[3].c[0].real; ai = a->d[3].c[0].imag;
  sum += ar * ar + ai * ai;
  ar = a->d[3].c[1].real; ai = a->d[3].c[1].imag;
  sum += ar * ar + ai * ai;
  ar = a->d[3].c[2].real; ai = a->d[3].c[2].imag;
  sum += ar * ar + ai * ai;
#endif
  return sum;
}

void magsq_wvec_sum(wilson_vector *a, Real *c) {
#ifndef FAST
  register int i;
  for (i = 0; i < 4; i++)
    magsq_su3vec_sum(&(a->d[i]), c);

#else  // FAST version for NCOL = DIMF = 3
  register Real ar, ai;

  ar = a->d[0].c[0].real; ai = a->d[0].c[0].imag;
  *c += ar * ar + ai * ai;
  ar = a->d[0].c[1].real; ai = a->d[0].c[1].imag;
  *c += ar * ar + ai * ai;
  ar = a->d[0].c[2].real; ai = a->d[0].c[2].imag;
  *c += ar * ar + ai * ai;

  ar = a->d[1].c[0].real; ai = a->d[1].c[0].imag;
  *c += ar * ar + ai * ai;
  ar = a->d[1].c[1].real; ai = a->d[1].c[1].imag;
  *c += ar * ar + ai * ai;
  ar = a->d[1].c[2].real; ai = a->d[1].c[2].imag;
  *c += ar * ar + ai * ai;

  ar = a->d[2].c[0].real; ai = a->d[2].c[0].imag;
  *c += ar * ar + ai * ai;
  ar = a->d[2].c[1].real; ai = a->d[2].c[1].imag;
  *c += ar * ar + ai * ai;
  ar = a->d[2].c[2].real; ai = a->d[2].c[2].imag;
  *c += ar * ar + ai * ai;

  ar = a->d[3].c[0].real; ai = a->d[3].c[0].imag;
  *c += ar * ar + ai * ai;
  ar = a->d[3].c[1].real; ai = a->d[3].c[1].imag;
  *c += ar * ar + ai * ai;
  ar = a->d[3].c[2].real; ai = a->d[3].c[2].imag;
  *c += ar * ar + ai * ai;
#endif
}
// -----------------------------------------------------------------
