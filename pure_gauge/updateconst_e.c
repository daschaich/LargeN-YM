// -----------------------------------------------------------------
// Update lattice with over-relaxed quasi-heat bath (qhb)
#include "pg_includes.h"

void updateconst_e(double Eint, double a) {
  // Check unitarity before doing anything
  check_unitarity();

#ifdef HMCLLR
  // Do HMC updates if specified
  update_hmc_const(Eint, a);
#else
  // Otherwise do over-relaxation and quasi-heat bath steps
  relax();
  monteconst_e(Eint, a);
#endif
  // Reunitarize the gauge field
  reunitarize();
}
// -----------------------------------------------------------------
