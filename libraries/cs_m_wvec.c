// -----------------------------------------------------------------
// Complex scalar multiplication on wilson_vector
// c <-- c + s * b
// c <-- a + s * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void c_scalar_mult_sum_wvec(wilson_vector *b, complex *s, wilson_vector *c) {
  register int i;
  for (i = 0; i < 4; i++)
    c_scalar_mult_sum_su3vec(&(b->d[i]), s, &(c->d[i]));
}

void c_scalar_mult_add_wvec(wilson_vector *a, wilson_vector *b,
                            complex *s, wilson_vector *c) {

  register int i;
  for (i = 0; i < 4; i++)
    c_scalar_mult_add_su3vec(&(a->d[i]), &(b->d[i]), s, &(c->d[i]));
}
// -----------------------------------------------------------------
