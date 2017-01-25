// -----------------------------------------------------------------
// Scalar multiplication on wilson_vector
// c <-- s * b
// c <-- c + s * b
// c <-- c - s * b
// c <-- a + s * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void scalar_mult_wvec(wilson_vector *b, Real s, wilson_vector *c) {
  register int i;
  for (i = 0; i < 4; i++)
    scalar_mult_vector(&(b->d[i]), s, &(c->d[i]));
}

void scalar_mult_sum_wvec(wilson_vector *b, Real s, wilson_vector *c) {
  register int i;
  for (i = 0; i < 4; i++)
    scalar_mult_sum_vector(&(b->d[i]), s, &(c->d[i]));
}

void scalar_mult_dif_wvec(wilson_vector *b, Real s, wilson_vector *c) {
  register int i;
  for (i = 0; i < 4; i++)
    scalar_mult_dif_vector(&(b->d[i]), s, &(c->d[i]));
}

void scalar_mult_add_wvec(wilson_vector *a, wilson_vector *b, Real s,
                          wilson_vector *c) {

  register int i;
  for (i = 0; i < 4; i++)
    scalar_mult_add_vector(&(a->d[i]), &(b->d[i]), s, &(c->d[i]));
}
// -----------------------------------------------------------------
