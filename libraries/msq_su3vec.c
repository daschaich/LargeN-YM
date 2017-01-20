// -----------------------------------------------------------------
// Squared magnitude of irrep vector
// Return ReTr[adag.a]
// c <-- c + ReTr[adag.a]
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

Real magsq_su3vec(su3_vector *a) {
  register Real sum = 0.0;
#ifndef FAST
  register int i;
  for (i = 0; i < DIMF; i++)
    sum += a->c[i].real * a->c[i].real + a->c[i].imag * a->c[i].imag;

#else  // FAST version for NCOL = DIMF = 3
  register Real temp;
  temp = a->c[0].real*a->c[0].real; sum += temp;
  temp = a->c[0].imag*a->c[0].imag; sum += temp;
  temp = a->c[1].real*a->c[1].real; sum += temp;
  temp = a->c[1].imag*a->c[1].imag; sum += temp;
  temp = a->c[2].real*a->c[2].real; sum += temp;
  temp = a->c[2].imag*a->c[2].imag; sum += temp;
#endif
  return sum;
}

void magsq_su3vec_sum(su3_vector *a, Real *c) {
#ifndef FAST
  register int i;
  for (i = 0; i < DIMF; i++)
    *c += a->c[i].real * a->c[i].real + a->c[i].imag * a->c[i].imag;

#else  // FAST version for NCOL = DIMF = 3
  register Real temp;
  temp = a->c[0].real*a->c[0].real; *c += temp;
  temp = a->c[0].imag*a->c[0].imag; *c += temp;
  temp = a->c[1].real*a->c[1].real; *c += temp;
  temp = a->c[1].imag*a->c[1].imag; *c += temp;
  temp = a->c[2].real*a->c[2].real; *c += temp;
  temp = a->c[2].imag*a->c[2].imag; *c += temp;
#endif
}
// -----------------------------------------------------------------
