// -----------------------------------------------------------------
// Adjoint of matrix
// b <-- adag
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void adjoint( matrix *a, matrix *b ){
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++)
      CONJG(a->e[j][i], b->e[i][j]);
  }
}
// -----------------------------------------------------------------
