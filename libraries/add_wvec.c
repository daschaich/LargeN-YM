// -----------------------------------------------------------------
// Add two wilson_vectors
// c <-- c + b
// c <-- a + b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void sum_wvec(wilson_vector *b, wilson_vector *c) {
   register int i;
   for (i = 0; i < 4; i++)
     sum_su3_vector(&(b->d[i]), &(c->d[i]));
}

void add_wilson_vector(wilson_vector *a, wilson_vector *b, wilson_vector *c) {
   register int i;
   for (i = 0; i < 4; i++)
     add_su3_vector(&(a->d[i]), &(b->d[i]), &(c->d[i]));
}
// -----------------------------------------------------------------
