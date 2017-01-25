// -----------------------------------------------------------------
// Clear the given irrep matrix
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void clear_mat(matrix *m) {
  register int i, j;
  for (i = 0; i < DIMF; i++) {
    for (j = 0; j < DIMF; j++) {
      m->e[i][j].real = 0.0;
      m->e[i][j].imag = 0.0;
    }
  }
}
// -----------------------------------------------------------------
