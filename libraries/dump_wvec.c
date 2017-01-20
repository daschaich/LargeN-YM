// -----------------------------------------------------------------
// Print the given wilson_vector
#include "../include/config.h"
#include <stdio.h>
#include "../include/complex.h"
#include "../include/su3.h"

void dump_wvec(wilson_vector *v) {
  register int i, j;
  for (i = 0; i < 4; i++) {
    for (j = 0; j < DIMF; j++)
      printf("  (%.4g, %.4g)", v->d[i].c[j].real, v->d[i].c[j].imag);
    printf("\n");
  }
}
// -----------------------------------------------------------------
