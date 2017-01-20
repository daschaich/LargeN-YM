// -----------------------------------------------------------------
// Return complex trace of the given fundamental matrix
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

/* Complex trace of an SU3 matrix */
complex trace_su3_f( su3_matrix_f *a ) {
  register complex tc;
  CADD(a->e[0][0], a->e[1][1], tc);
#if NCOL > 2
  CSUM(tc, a->e[2][2]);
#if NCOL > 3
  CSUM(tc, a->e[3][3]);
#endif
#endif
  return tc;
}
// -----------------------------------------------------------------
