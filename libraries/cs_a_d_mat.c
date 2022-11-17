// -----------------------------------------------------------------
// Add complex scalar to each diagonal element of matrix
// a <-- a + c * I
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void c_scalar_add_diag(matrix *a, complex *c) {
  register int i;
  for (i = 0; i < NCOL; i++) {
    a->e[i][i].real += c->real;
    a->e[i][i].imag += c->imag;
  }
}
