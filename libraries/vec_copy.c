// -----------------------------------------------------------------
// Copy an irrep vector
// b <-- a
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void vec_copy(vector *a, vector *b) {
  register int i;
  for (i = 0; i < DIMF; i++) {
    b->c[i].real = a->c[i].real;
    b->c[i].imag = a->c[i].imag;
  }
}
// -----------------------------------------------------------------
