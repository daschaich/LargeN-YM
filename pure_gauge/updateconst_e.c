// -----------------------------------------------------------------
// Update lattice with microcanonical over-relaxed quasi-heat bath (qhb)
#include "pg_includes.h"

int updateconst_e(double Eint, double a) {
  int iters=0;

  // Check unitarity before doing anything
  check_unitarity();

  // Do over-relaxation and quasi-heat bath steps
  relax(steps);
  monteconst_e(stepsQ, Eint, a);

  // Reunitarize the gauge field
  reunitarize();

  if (steps > 0)
    return (iters / steps);
  else
    return(-99);
}
// -----------------------------------------------------------------
