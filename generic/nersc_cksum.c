/* Two utilities used in the NERSC archive formate */

/* Compute the low order 32 bits of the unsigned integer sum of the
   float precision real and complex parts of the elements of the gauge
   matrices.
*/

/* Computes the mean global sum of the trace of the gauge links --
   used to aid checking lattice file integrity */

#include "generic_includes.h"

u_int32type nersc_cksum() {
  return 0;
}

void linktrsum(double_complex *linktr) {
  int i, dir;
  site *s;
  matrix_f *a;

  *linktr = dcmplx(0.0, 0.0);
  FORALLSITES(i, s) {
    FORALLUPDIR(dir) {
      a = &s->linkf[dir];
      CSUM(*linktr, a->e[0][0]);
      CSUM(*linktr, a->e[1][1]);
#if NCOL > 2
      CSUM(*linktr, a->e[2][2]);
#if NCOL > 3
      CSUM(*linktr, a->e[3][3]);
#endif
#endif
    }
  }

  g_dcomplexsum(linktr);
  CDIVREAL(*linktr, (4 * volume), *linktr);
}
