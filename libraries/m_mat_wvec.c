// -----------------------------------------------------------------
// Irrep matrix multiplying Wilson vector
// c <-- a.b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void mult_mat_wilson_vec(su3_matrix *a, wilson_vector *b,
                         wilson_vector *c) {

    register int i;
    for (i = 0; i < 4; i++)
      mult_su3_mat_vec(a, &(b->d[i]), &(c->d[i]));
}
// -----------------------------------------------------------------
