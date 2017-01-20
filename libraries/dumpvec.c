// -----------------------------------------------------------------
// Print the given irrep vector
#include "../include/config.h"
#include <stdio.h>
#include "../include/complex.h"
#include "../include/su3.h"

void dump_vec(su3_vector *v) {
  register int j;
  for (j = 0; j < DIMF; j++)
    printf("  (%.4g, %.4g)", v->c[j].real, v->c[j].imag);
}
// -----------------------------------------------------------------
