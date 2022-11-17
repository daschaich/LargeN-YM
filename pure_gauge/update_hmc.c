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
int update_hmc() {
  int step, iters = 0;
  Real xrandom, tr;
  Real eps = traj_length / (Real)hmc_steps;
  double fnorm = 0.0, startaction = 0.0, endaction, change;

  // Refresh the momenta
  ranmom();

  // Find initial action
  startaction = action();

  // Copy link field to old_link
  gauge_field_copy(F_OFFSET(link[0]), F_OFFSET(old_link[0]));

  // Do microcanonical updating
  // First u(t/2)
  update_u(0.5 * eps);

  // Inner steps p(t) u(t)
  for (step = 0; step < hmc_steps; step++) {
    tr = update_h(eps);
    fnorm += tr;
    if (tr > max_f)
      max_f = tr;

    if (step < hmc_steps - 1)
      update_u(eps);
    else
      update_u(0.5 * eps);    // Final u(t/2)
  }

  // Reunitarize the gauge field
  reunitarize();

  // Find ending action
  endaction = action();
  change = endaction - startaction;

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
  if (exp(-change) < (double)xrandom) {
    if (traj_length > 0.0)
      gauge_field_copy(F_OFFSET(old_link[0]), F_OFFSET(link[0]));

    node0_printf("REJECT: delta S = %.4g start S = %.12g end S = %.12g\n",
                 change, startaction, endaction);
  }
  else {
    node0_printf("ACCEPT: delta S = %.4g start S = %.12g end S = %.12g\n",
                 change, startaction, endaction);
  }

  if (traj_length > 0.0) {
    node0_printf("IT_PER_TRAJ %d\n", iters);
    node0_printf("MONITOR_FORCE %.4g %.4g\n",
                 fnorm / (double)(2 * hmc_steps), max_f);
    return iters;
  }
  else
    return(-99);
}

// TODO: Should be able to use update_hmc with modified action()...
int update_hmc_const(double Eint, double a) {
#ifdef LLR
#if 0
  int step, iters = 0;
  Real xrandom, tr;
  Real eps = traj_length / (Real)hmc_steps;
  double fnorm = 0.0, startactionHMC = 0.0, endactionHMC, change;
  double startlatticeaction = 0.0;
  double endlatticeaction;
  double regulardelta;
  // Refresh the momenta
  ranmom();

  // Find initial action
  startactionHMC = action_HMC(a);
  startlatticeaction = U_action();
  // Copy link field to old_link
  gauge_field_copy(F_OFFSET(link[0]), F_OFFSET(old_link[0]));

  // Do microcanonical updating
  // First u(t/2)
  update_u(0.5 * eps);

  // Inner steps p(t) u(t)
  for (step = 0; step < hmc_steps; step++) {
    tr = update_h_const(eps,Eint,a);
    fnorm += tr;
    if (tr > max_f)
      max_f = tr;

    if (step < hmc_steps - 1)
      update_u(eps);
    else
      update_u(0.5 * eps);    // Final u(t/2)
  }

  // Reunitarize the gauge field
  reunitarize();

  // Find ending action
  endactionHMC = action_HMC(a);
  endlatticeaction = U_action();
  change = endactionHMC - startactionHMC;
  regulardelta = endlatticeaction - startlatticeaction;

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
  //if (exp(-change) < (double)xrandom) {
  if ((double)xrandom < exp(-change-pow(regulardelta,2.0)/(2.0*pow(delta,2.0))-regulardelta*(startlatticeaction-Eint-delta*0.5)/pow(delta,2.0))) {
    if (traj_length > 0.0)
      gauge_field_copy(F_OFFSET(old_link[0]), F_OFFSET(link[0]));

    node0_printf("REJECT: delta S = %.4g start S = %.12g end S = %.12g\n",
                 change, startactionHMC, endactionHMC);
  }
  else {
    node0_printf("ACCEPT: delta S = %.4g start S = %.12g end S = %.12g\n",
                 change, startactionHMC, endactionHMC);
  }

  if (traj_length > 0.0) {
    node0_printf("IT_PER_TRAJ %d\n", iters);
    node0_printf("MONITOR_FORCE %.4g %.4g\n",
                 fnorm / (double)(2 * hmc_steps), max_f);
    return iters;
  }
  else
#endif
#endif
    return(-99);
}
// -----------------------------------------------------------------
