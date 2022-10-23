// -----------------------------------------------------------------
// Copy a matrix
// b <-- a
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void mat_copy(matrix *a, matrix *b) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      b->e[i][j].real = a->e[i][j].real;
      b->e[i][j].imag = a->e[i][j].imag;
    }
  }
}
// -----------------------------------------------------------------
