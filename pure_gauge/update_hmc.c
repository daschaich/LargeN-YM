// -----------------------------------------------------------------
// Update lattice with simple leapfrog MD integrator
// Begin at "integral" time, with H and U evaluated at the same time
#include "pg_includes.h"
#ifdef HAVE_IEEEFP_H
#include <ieeefp.h>         // For "finite"
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Helper function copies a gauge field as an array of four matrices
void gauge_field_copy(field_offset src, field_offset dest) {
  register int i, dir, src2, dest2;
  register site *s;

  FORALLSITES(i, s) {
    src2 = src;
    dest2 = dest;
    FORALLUPDIR(dir) {
      mat_copy((matrix *)F_PT(s, src2), (matrix *)F_PT(s, dest2));
      src2 += sizeof(matrix);
      dest2 += sizeof(matrix);
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void update_hmc(double E_min) {
  int step;
  Real xrandom, tr;
  Real eps = traj_length / (Real)hmc_steps;
  double fnorm = 0.0, startaction = 0.0, endaction, change;
  
  double deltaH2 = 0;
  double regularstart = 0;
  double regularend = 0;
  double regulardelta = 0;
  
  // Refresh the momenta
  ranmom();

  // Find initial action
  startaction = action(E_min);
  regularstart = gauge_action();
  // Copy link field to old_link
  gauge_field_copy(F_OFFSET(link[0]), F_OFFSET(old_link[0]));

  // Do microcanonical updating
  // First u(t/2)
  update_u(0.5 * eps);

  // Inner steps p(t) u(t)
  for (step = 0; step < hmc_steps; step++) {
    tr = update_h(eps,E_min);
    fnorm += tr;
    if (tr > max_f)
      max_f = tr;

    if (step < hmc_steps - 1)
      update_u(eps);
    else
      update_u(0.5 * eps);    // Final u(t/2)
    
    node0_printf("gauge_action after update= %.4g\n", gauge_action());
  }

  // Reunitarize the gauge field
  reunitarize();

  // Find ending action
  endaction = action(E_min);
  regularend = gauge_action();
  change = endaction - startaction;
  regulardelta = regularend - regularstart;
  // Reject configurations giving overflow
#ifndef HAVE_IEEEFP_H
  if (fabs((double)change) > 1e20) {
#else
  if (!finite((double)change)) {
#endif
    node0_printf("WARNING: Correcting Apparent Overflow: Delta S = %.4g\n",
                 change);
    change = 1.0e20;
  }

  // Decide whether to accept, if not, copy old link field back
  // Careful -- must generate only one random number for whole lattice
  if (this_node == 0)
    xrandom = myrand(&node_prn);
  broadcast_float(&xrandom);

if (constrained == 0) {
  if (exp(-change) < (double)xrandom) {
    if (traj_length > 0.0)
      gauge_field_copy(F_OFFSET(old_link[0]), F_OFFSET(link[0]));

    node0_printf("REJECT: delta S = %.4g start S = %.12g end S = %.12g\n",
                 change, startaction, endaction);                 
                 deltaH2 = deltaH2 + pow(-change,2); //-pow(regulardelta,2.0)/(2.0*pow(delta,2.0))-regulardelta*(regularstart-E_min-delta*0.5)/pow(delta,2.0),2);
  }
  else {
    node0_printf("ACCEPT: delta S = %.4g start S = %.12g end S = %.12g\n",
                 change, startaction, endaction);
                 deltaH2 = deltaH2 + pow(-change,2); //-pow(regulardelta,2.0)/(2.0*pow(delta,2.0))-regulardelta*(regularstart-E_min-delta*0.5)/pow(delta,2.0),2);
  }
}

else {
 if ((double)xrandom < exp(-change-pow(regulardelta,2.0)/(2.0*pow(delta,2.0))-regulardelta*(regularstart-E_min-delta*0.5)/pow(delta,2.0))) {
  if (traj_length > 0.0)
      node0_printf("ACCEPT: delta S = %.4g start S = %.12g end S = %.12g\n",
                 change, startaction, endaction);
                 deltaH2 = deltaH2 + pow(-change-pow(regulardelta,2.0)/(2.0*pow(delta,2.0))-regulardelta*(regularstart-E_min-delta*0.5)/pow(delta,2.0),2);
 }
 else {
    if (traj_length > 0.0)
      gauge_field_copy(F_OFFSET(old_link[0]), F_OFFSET(link[0]));

    node0_printf("REJECT: delta S = %.4g start S = %.12g end S = %.12g\n",
                 change, startaction, endaction);                 
                 deltaH2 = deltaH2 + pow(-change-pow(regulardelta,2.0)/(2.0*pow(delta,2.0))-regulardelta*(regularstart-E_min-delta*0.5)/pow(delta,2.0),2);
 }
}
#ifdef LLR
      if (constrained == 1) {
        node0_printf("DeltaLLR_H^2 = %.12g\n", deltaH2);
      }
      else {
        node0_printf("Delta_H^2 = %.12g\n", deltaH2);
      }
#endif
  if (traj_length > 0.0) {
    node0_printf("MONITOR_FORCE %.4g %.4g\n",
                 fnorm / (double)(2 * hmc_steps), max_f);
  }
}
// -----------------------------------------------------------------
