// -----------------------------------------------------------------
// Update lattice, choose between hybrid Monte Carlo trajectories
// or over-relaxed quasi-heat bath (qhb) sweeps
#include "pg_includes.h"

void update() {
  // Check unitarity before doing anything
  check_unitarity();

#ifdef HMC_ALGORITHM
  // Do HMC updates if specified
  
#else
  // Otherwise do over-relaxation and quasi-heat bath steps
  relax();
  monte();
#endif

  // Reunitarize the gauge field
  reunitarize();
}
// -----------------------------------------------------------------
