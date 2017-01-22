// -----------------------------------------------------------------
// Clear a fundamental vector
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void clearvec_f(su3_vector_f *v) {
  v->c[0].real = 0.0;
  v->c[0].imag = 0.0;
  v->c[1].real = 0.0;
  v->c[1].imag = 0.0;
#if NCOL > 2
  v->c[2].real = 0.0;
  v->c[2].imag = 0.0;
#if NCOL > 3
  v->c[3].real = 0.0;
  v->c[3].imag = 0.0;
#if NCOL > 4
  int i;
  for (i = 4; i < NCOL; i++) {
    v->c[i].real = 0.0;
    v->c[i].imag = 0.0;
  }
#endif
#endif
#endif
}
// -----------------------------------------------------------------
