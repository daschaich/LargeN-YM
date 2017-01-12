// -----------------------------------------------------------------
// Return real part of dot product of two irrep vectors
// ReTr[adag.b]
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

Real su3_rdot(su3_vector *a, su3_vector *b) {
  register Real tr,sum;
#if DIMF > 6
  register int i;
#endif

  sum = a->c[0].real * b->c[0].real;
  tr = a->c[0].imag * b->c[0].imag; sum += tr;
  tr = a->c[1].real * b->c[1].real; sum += tr;
  tr = a->c[1].imag * b->c[1].imag; sum += tr;
#if DIMF > 2
  tr = a->c[2].real * b->c[2].real; sum += tr;
  tr = a->c[2].imag * b->c[2].imag; sum += tr;
#if DIMF > 3
  tr = a->c[3].real * b->c[3].real; sum += tr;
  tr = a->c[3].imag * b->c[3].imag; sum += tr;
#if DIMF > 4
  tr = a->c[4].real * b->c[4].real; sum += tr;
  tr = a->c[4].imag * b->c[4].imag; sum += tr;
#if DIMF > 5
  tr = a->c[5].real * b->c[5].real; sum += tr;
  tr = a->c[5].imag * b->c[5].imag; sum += tr;
#if DIMF > 6
  for (i = 6; i < DIMF; i++) {
    tr = a->c[i].real * b->c[i].real; sum += tr;
    tr = a->c[i].imag * b->c[i].imag; sum += tr;
  }
#endif
#endif
#endif
#endif
#endif

  return sum;
}
