// -----------------------------------------------------------------
// Wrapper for Wilson and Wilson-clover quark inverter
// The source routine must be called before this inversion routine
#include "generic_wilson_includes.h"

// src, dest and sav are all wilson_vectors
// dest contains first the initial guess and then the answer
// sav is for saving the source
// qic is the quark inverter control structure
// dmp allows us to pass Dirac matrix parameters
// Return value is number of iterations
int wilson_invert(field_offset src, field_offset dest, field_offset sav,
                  int (*invert_func)(field_offset src, field_offset dest,
                                     quark_invert_control *qic, void *dmp),
                  quark_invert_control *qic, void *dmp) {

  register int i;
  register site *s;
  int tot_iters, irestart, Minsav = qic->min;

  // Store the source for future use
  FORALLSITES(i, s)
    copy_wvec((wilson_vector *)F_PT(s, src), (wilson_vector *)F_PT(s, sav));

  // Inversion with restart (appropriate for CG and BiCG)
  for (tot_iters = 0, irestart = 0; irestart < qic->nrestart; irestart++) {
    tot_iters += invert_func(src, dest, qic, dmp);

    // Check for convergence
    if (qic->size_r < qic->resid)
      break;

    // Restore the source for restart (but not for exit)
    FORALLSITES(i, s)
      copy_wvec((wilson_vector *)F_PT(s, sav), (wilson_vector *)F_PT(s, src));

    // Restart trial solution is not zero any more
    qic->start_flag = 1;

    // No minimum when restarting
    qic->min = 0;
  }

  // Warn if inversion fails
  if (qic->size_r > qic->resid) {
    node0_printf("wilson_invert did not converge after %d iterations, ",
                 tot_iters);
    node0_printf("size_r = %.8g wanted %.8g\n", qic->size_r, qic->resid);
  }

  // Restore minimum
  qic->min = Minsav;
  return tot_iters;
}
// -----------------------------------------------------------------
