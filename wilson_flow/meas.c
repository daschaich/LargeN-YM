// -----------------------------------------------------------------
// Measure observables after flowing for time t
// Prepend flow time to output
#include "wflow_includes.h"

void meas(Real t) {
  int j, k;
  double norm, rr[10];
  complex plp, xplp;

  // Wilson loops
  gauge_loops(rr);
  for (j = 0; j < nloop; j++) {
    norm = 1.0 / (double)(volume * loop_num[j]);
    for (k = 0; k < nreps; k++) {
      rr[j + k * nloop] *= norm;
      node0_printf("LOOPS %g %d %d 0 0 %.8g\n", t, j, k, rr[j + k * nloop]);
    }
  }

  // Polyakov loop and one spatial Wilson line
  // Keep "ORIG" in output for analysis script backward compatibility
  // First call will print PL(x, y, z) if LOCALPOLY defined
  plp = ploop(TUP);
  xplp = ploop(XUP);
  node0_printf("POLYA ORIG %g %.6g %.6g %.6g %.6g\n",
               t, plp.real, plp.imag, xplp.real, xplp.imag);
}
// -----------------------------------------------------------------
