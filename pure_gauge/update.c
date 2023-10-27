// -----------------------------------------------------------------
// Update lattice, choosing between hybrid Monte Carlo trajectories
// or over-relaxed quasi-heatbath sweeps
// In both cases Check unitarity before doing anything
// and reunitarize gauge field at end
#include "pg_includes.h"

#if defined(ORA) || defined(LLR)
void update_ora() {
  check_unitarity();
  relax();
  monte();
  reunitarize();
}
#endif

// HMC updates including energy interval info for LLR
#ifdef HMC
void update_hmc(Real C, double E_min) {
  check_unitarity();
  hmc_traj(C, E_min);
  reunitarize();
}
#endif
// -----------------------------------------------------------------
