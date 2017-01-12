// -----------------------------------------------------------------
// Scalar multiplication on irrep vector
// b <- s * a
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void scalar_mult_su3_vector(su3_vector *a, Real s, su3_vector *b) {
  register int i;
  for (i = 0; i < DIMF; i++) {
    b->c[i].real = s * a->c[i].real;
    b->c[i].imag = s * a->c[i].imag;
  }
}
// -----------------------------------------------------------------
