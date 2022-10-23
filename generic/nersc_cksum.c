// -----------------------------------------------------------------
// Two utilities used to check saved lattice files
#include "generic_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Turn this one off for now
u_int32type nersc_cksum() {
  return 0;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute the mean global sum of the trace of the gauge links
// Use to check lattice file integrity
void linktrsum(double_complex *linktr) {
  int i, dir;
  site *s;
  matrix *a;

  *linktr = dcmplx(0.0, 0.0);
  FORALLSITES(i, s) {
    FORALLUPDIR(dir) {
      a = &(s->linkf[dir]);
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
  CDIVREAL(*linktr, (4.0 * volume), *linktr);
}
// -----------------------------------------------------------------
