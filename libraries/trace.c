// -----------------------------------------------------------------
// Complex trace of the given matrix
// Return Tr[m]
// c <-- c + Tr[m]
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

complex trace(matrix *m) {
  register complex tc;
  CADD(m->e[0][0], m->e[1][1], tc);
#if NCOL > 2
  CSUM(tc, m->e[2][2]);
#if NCOL > 3
  CSUM(tc, m->e[3][3]);
#if NCOL > 4
  register int i;
  for (i = 4; i < NCOL; i++)
    CSUM(tc, m->e[i][i]);
#endif
#endif
#endif
  return tc;
}

void trace_sum(matrix *m, complex *c) {
  CSUM(*c, m->e[0][0]);
  CSUM(*c, m->e[1][1]);
#if NCOL > 2
  CSUM(*c, m->e[2][2]);
#if NCOL > 3
  CSUM(*c, m->e[3][3]);
#if NCOL > 4
  register int i;
  for (i = 4; i < NCOL; i++)
    CSUM(*c, m->e[i][i]);
#endif
#endif
#endif
}
// -----------------------------------------------------------------
