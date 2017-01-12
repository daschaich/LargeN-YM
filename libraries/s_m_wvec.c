// -----------------------------------------------------------------
// Scalar multiplication on wilson_vector
// b <-- s * a
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void scalar_mult_wvec(wilson_vector *a, Real s, wilson_vector *b) {
  register int i;
  for (i = 0; i < 4; i++)
    scalar_mult_su3_vector(&(a->d[i]), s, &(b->d[i]));
}
// -----------------------------------------------------------------
