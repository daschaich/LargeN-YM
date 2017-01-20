// -----------------------------------------------------------------
// Print the given irrep matrix
#include "../include/config.h"
#include <stdio.h>
#include "../include/complex.h"
#include "../include/su3.h"

void dumpmat(su3_matrix *m) {
  int i, j;
  for (i = 0; i < DIMF; i++) {
    for (j = 0; j < DIMF; j++)
      printf("  (%.4g, %.4g)", m->e[i][j].real, m->e[i][j].imag);
    printf("\n");
  }
}
// -----------------------------------------------------------------
