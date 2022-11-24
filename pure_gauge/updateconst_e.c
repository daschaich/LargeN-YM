// -----------------------------------------------------------------
// Update lattice
// Choose between hybrid Monte Carlo trajectories
// or over-relaxed quasi-heatbath (qhb) sweeps
#include "pg_includes.h"

void updateconst_e(double Eint, double a) {
#ifdef LLR
  // Check unitarity before doing anything
  check_unitarity();

#ifdef HMC
  // Do HMC updates if specified
  update_hmc();
#else
  // Otherwise do over-relaxation and quasi-heatbath steps
  relax();
  monteconst_e();
#endif

  // Reunitarize the gauge field
  reunitarize();
#endif
}
// -----------------------------------------------------------------
