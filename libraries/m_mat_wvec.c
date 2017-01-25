// -----------------------------------------------------------------
// Irrep matrix multiplying Wilson vector
// c <-- a.b
// c <-- adag.b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void mult_mat_wvec(matrix *a, wilson_vector *b, wilson_vector *c) {
  register int i;
  for (i = 0; i < 4; i++)
    mult_mat_vec(a, &(b->d[i]), &(c->d[i]));
}

void mult_adj_mat_wvec(matrix *a, wilson_vector *b, wilson_vector *c) {
  register int i;
  for (i = 0; i < 4; i++)
    mult_adj_mat_vec(a, &(b->d[i]), &(c->d[i]));
}
// -----------------------------------------------------------------
