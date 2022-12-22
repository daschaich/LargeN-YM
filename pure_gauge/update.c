// -----------------------------------------------------------------
// Update lattice, choosing between hybrid Monte Carlo trajectories
// or over-relaxed quasi-heatbath (qhb) sweeps
#include "pg_includes.h"

void update() {
  // Check unitarity before doing anything
  check_unitarity();

#ifdef HMC
  // Do HMC updates if specified
  update_hmc(0.0);
#else
  // Otherwise do over-relaxation and quasi-heatbath steps
  relax();
  monte();
#endif

  // Reunitarize the gauge field
  reunitarize();
}
// -----------------------------------------------------------------
