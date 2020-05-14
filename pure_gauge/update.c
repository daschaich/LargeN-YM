// -----------------------------------------------------------------
// Update lattice with microcanonical over-relaxed quasi-heat bath (qhb)
#include "pg_includes.h"

int update() {
  int iters=0;

  // Check unitarity before doing anything
  check_unitarity();

  // Do over-relaxation and quasi-heat bath steps
  relax(steps);
  monte(stepsQ);

  // Reunitarize the gauge field
  reunitarize();

  if (steps > 0)
    return (iters / steps);
  else
    return(-99);
}
// -----------------------------------------------------------------
