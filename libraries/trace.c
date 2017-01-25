// -----------------------------------------------------------------
// Return complex trace of the given irrep matrix
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

complex trace(matrix *a) {
  register complex tc;
  CADD(a->e[0][0], a->e[1][1], tc);
#if DIMF > 2
  CSUM(tc, a->e[2][2]);
#if DIMF > 3
  CSUM(tc, a->e[3][3]);
#if DIMF > 4
  CSUM(tc, a->e[4][4]);
#if DIMF > 5
  CSUM(tc, a->e[5][5]);
#if DIMF > 6
  register int i;
  for (i = 6; i < DIMF; i++)
    CSUM(tc, a->e[i][i]);
#endif
#endif
#endif
#endif
#endif
  return tc;
}
// -----------------------------------------------------------------
