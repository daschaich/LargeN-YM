// -----------------------------------------------------------------
// Clear the given vector
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void clearvec(vector *v) {
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
  v->c[4].real = 0.0;
  v->c[4].imag = 0.0;
#if NCOL > 5
  v->c[5].real = 0.0;
  v->c[5].imag = 0.0;
#if NCOL > 6
  register int i;
  for (i = 6; i < NCOL; i++) {
    v->c[i].real = 0.0;
    v->c[i].imag = 0.0;
  }
#endif
#endif
#endif
#endif
#endif
}
// -----------------------------------------------------------------
