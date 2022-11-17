// -----------------------------------------------------------------
// Add scalar to each diagonal element of matrix
// a <-- a + s * I
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void scalar_add_diag(matrix *a, Real s) {
  register int i;
  for (i = 0; i < NCOL; i++)
    a->e[i][i].real += s;
}
// -----------------------------------------------------------------
