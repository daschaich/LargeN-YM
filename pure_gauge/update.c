// -----------------------------------------------------------------
// Update lattice with over-relaxed quasi-heat bath (qhb)
#include "pg_includes.h"

void update() {
  // Check unitarity before doing anything
  check_unitarity();

  // Do over-relaxation and quasi-heat bath steps
  relax();
  monte();

  // Reunitarize the gauge field
  reunitarize();
}
// -----------------------------------------------------------------
