// -----------------------------------------------------------------
// Update lattice with over-relaxed quasi-heat bath (qhb)
#include "pg_includes.h"

void updateconst_e(double Eint, double a) {
  // Check unitarity before doing anything
  check_unitarity();

  // Do over-relaxation and quasi-heat bath steps
  relax();
  monteconst_e(Eint, a);

  // Reunitarize the gauge field
  reunitarize();
}
// -----------------------------------------------------------------
