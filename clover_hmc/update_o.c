// -----------------------------------------------------------------
// Update lattice, with Omelyan integrator and a Hasenbusch mass

// Note that for the final accept/reject,
// we already have a good solution to the CG:
// The last update was of the momenta

// The input file gave us the number of levels, trajectory length, and
// nsteps[MAX_MASSES + 1] (the gauge level is the last)
// num_masses <= MAX_MASSES is the total number of
// real dynamical fermion doublets plus the Hasenbusch preconditioners

// This routine begins at "integral" time
// with H and U evaluated at same time
#include "cl_dyn_includes.h"

// Local function prototypes
void predict_next_psi(Real *oldtime, Real *newtime,
                      Real *nexttime, int level);
int update_step(Real *oldtime, Real *newtime, Real *nexttime,
                double *fnorm, double *gnorm);
double returntrlogA;

#ifdef LU
#define FORMYSITESDOMAIN FOREVENSITESDOMAIN
#else
#define FORMYSITESDOMAIN FORALLSITESDOMAIN
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Omelyan gauge updates -- `dirty' sped-up version
// Return norm of force for averaging
double update_gauge_step(Real eps) {
  int n = nsteps[MAX_MASSES], i;
  double norm = 0.0;

#ifdef CG_DEBUG
  node0_printf("UPDATE: %d gauge steps of %.4g\n", n, eps);
#endif

  norm += gauge_force(eps * LAMBDA);
  for (i = 1; i <= n; i++) {
    update_u(0.5 * eps);
    norm += gauge_force(eps * LAMBDA_MID);
    update_u(0.5 * eps);

    // 2lambda step except on last iteration; then only lambda
    if (i < n)
      norm += gauge_force(eps * TWO_LAMBDA);
    else
      norm += gauge_force(eps * LAMBDA);
  }

  // Reunitarize the gauge field and update DIMFxDIMF matrices
  reunitarize();
  fermion_rep();
  return norm / n;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Use linear extrapolation to predict next conjugate gradient solution
// Only need even sites
void predict_next_psi(Real *oldtime, Real *newtime, Real *nexttime,
                      int level) {

  register int i;
  register site *s;
  register Real x;
  wilson_vector tvec;

  if (newtime[level] != oldtime[level]) {
    x = (nexttime[level] - newtime[level])
      / (newtime[level] - oldtime[level]);
  }
  else
    x = 0.0;

  if (oldtime[level] < 0.0) {
    FOREVENSITES(i, s)
      s->old_psi[level] = s->psi[level];
  }
  else  {
    FOREVENSITES(i,s) {
      sub_wilson_vector(&(s->psi[level]), &(s->old_psi[level]), &tvec);
      s->old_psi[level] = s->psi[level];
      scalar_mult_add_wvec(&(s->psi[level]), &tvec,x, &(s->psi[level]));
    }
  }
  oldtime[level] = newtime[level];
  newtime[level] = nexttime[level];
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Three-level Omelyan integrator for Hasenbusch-preconditioned code
int update_step(Real *old_cg_time, Real *cg_time, Real *next_cg_time,
                double *fnorm, double *gnorm)  {

  int iters = 0, outer, inner, level;
  Real CKU0 = kappa * clov_c / (u0 * u0 * u0), final_rsq;
  Real f_eps0, f_eps1, g_eps, mshift, tr;

  f_eps0 = traj_length / (Real)nsteps[0];
  f_eps1 = f_eps0 / (2.0 * (Real)nsteps[1]);
  g_eps = f_eps1 / (2.0 * (Real)nsteps[MAX_MASSES]);

  if (num_masses == 1) {
    level = 0;
    mshift = 0.0;
  }
  else {
    level = 1;
    mshift = shift;
  }
#ifdef CG_DEBUG
  node0_printf("level %d mshift %e\n", level, mshift);
#endif

  // Initial lambda update steps for both inner and outer force loops
  // (Already did CG on both levels to get starting action)
  tr = fermion_force(f_eps1 * LAMBDA, f_eps0 * LAMBDA);
  fnorm[0] += tr;
  if (tr > max_ff[0])
    max_ff[0] = tr;
  for (outer = 1; outer <= nsteps[0]; outer++) {
    for (inner = 1; inner <= nsteps[1]; inner++) {
      tr = update_gauge_step(g_eps);
      *gnorm += tr;
      if (tr > max_gf)
        max_gf = tr;

      // (1 - 2lambda) inner force update step
      next_cg_time[level] = cg_time[level] + f_eps1;
      predict_next_psi(old_cg_time, cg_time, next_cg_time, level);
      free_clov();
      make_clov(CKU0);
#ifdef LU
      returntrlogA = make_clovinv(ODD);
#endif
      iters += congrad_cl_m(niter, rsqmin, &final_rsq, F_OFFSET(chi[level]),
                            F_OFFSET(psi[level]), mshift);

      tr = fermion_force(f_eps1 * LAMBDA_MID, 0.0);
      fnorm[1] += tr;
      if (tr > max_ff[1])
        max_ff[1] = tr;

      tr = update_gauge_step(g_eps);
      *gnorm += tr;
      if (tr > max_gf)
        max_gf = tr;

      // 2lambda inner force update step, since combined
      // with initial update step of next loop
      if (inner < nsteps[1]) {
        next_cg_time[level] = cg_time[level] + f_eps1;
        predict_next_psi(old_cg_time, cg_time, next_cg_time, level);
        free_clov();
        make_clov(CKU0);
#ifdef LU
        returntrlogA = make_clovinv(ODD);
#endif
        iters += congrad_cl_m(niter, rsqmin, &final_rsq, F_OFFSET(chi[level]),
                              F_OFFSET(psi[level]), mshift);
        tr = fermion_force(f_eps1 * TWO_LAMBDA, 0.0);
        fnorm[1] += tr;
        if (tr > max_ff[1])
          max_ff[1] = tr;
      }
    } // Inner update steps

    // (1 - 2lambda) outer force update step
    next_cg_time[level] = cg_time[level] + f_eps1;
    predict_next_psi(old_cg_time, cg_time, next_cg_time,level);
    free_clov();
    make_clov(CKU0);
#ifdef LU
    returntrlogA = make_clovinv(ODD);
#endif
    iters += congrad_cl_m(niter, rsqmin, &final_rsq, F_OFFSET(chi[level]),
                          F_OFFSET(psi[level]), mshift);

    if (num_masses == 2) {
      next_cg_time[0] = cg_time[0] + f_eps0;
      predict_next_psi(old_cg_time, cg_time, next_cg_time, 0);
      iters += congrad_cl_m(niter, rsqmin, &final_rsq, F_OFFSET(chi[0]),
                            F_OFFSET(psi[0]), 0.0);
    }

    tr = fermion_force(f_eps1 * TWO_LAMBDA, f_eps0 * LAMBDA_MID);
    fnorm[0] += tr;
    if (tr > max_ff[0])
      max_ff[0] = tr;

    // Initial lambda update step done above
    for (inner = 1; inner <= nsteps[1]; inner++) {
      tr = update_gauge_step(g_eps);
      *gnorm += tr;
      if (tr > max_gf)
        max_gf = tr;

      // (1 - 2lambda) inner force update step
      next_cg_time[level] = cg_time[level] + f_eps1;
      predict_next_psi(old_cg_time, cg_time, next_cg_time, level);
      free_clov();
      make_clov(CKU0);
#ifdef LU
      returntrlogA = make_clovinv(ODD);
#endif
      iters += congrad_cl_m(niter, rsqmin, &final_rsq, F_OFFSET(chi[level]),
                            F_OFFSET(psi[level]), mshift);
      tr = fermion_force(f_eps1 * LAMBDA_MID, 0.0);
      fnorm[1] += tr;
      if (tr > max_ff[1])
        max_ff[1] = tr;

      tr = update_gauge_step(g_eps);
      *gnorm += tr;
      if (tr > max_gf)
        max_gf = tr;

      // 2lambda inner force update step, except on last iteration
      if (inner < nsteps[1]) {
        next_cg_time[level] = cg_time[level] + f_eps1;
        predict_next_psi(old_cg_time, cg_time, next_cg_time, level);
        free_clov();
        make_clov(CKU0);
#ifdef LU
        returntrlogA = make_clovinv(ODD);
#endif
        iters += congrad_cl_m(niter,rsqmin,&final_rsq,F_OFFSET(chi[level]),
            F_OFFSET(psi[level]),mshift);
        tr = fermion_force(f_eps1 * TWO_LAMBDA, 0.0);
        fnorm[1] += tr;
        if (tr > max_ff[1])
          max_ff[1] = tr;
      }
    } // Inner update steps

    // Except on last iteration, 2lambda update step for both
    // On last iteration, just lambda update step for both
    next_cg_time[level] = cg_time[level] + f_eps1;
    predict_next_psi(old_cg_time, cg_time, next_cg_time, level);
    free_clov();
    make_clov(CKU0);
#ifdef LU
    returntrlogA = make_clovinv(ODD);
#endif
    iters += congrad_cl_m(niter, rsqmin, &final_rsq, F_OFFSET(chi[level]),
                          F_OFFSET(psi[level]), mshift);
    if (num_masses == 2) {
      next_cg_time[0] = cg_time[0] + f_eps0;
      predict_next_psi(old_cg_time,cg_time,next_cg_time,0);
      iters += congrad_cl_m(niter, rsqmin, &final_rsq, F_OFFSET(chi[0]),
                            F_OFFSET(psi[0]), 0.0);
    }

    if (outer < nsteps[0])
      tr += fermion_force(f_eps1 * TWO_LAMBDA, f_eps0 * TWO_LAMBDA);
    else
      tr += fermion_force(f_eps1 * LAMBDA, f_eps0 * LAMBDA);
    fnorm[0] += tr;
    if (tr > max_ff[0])
      max_ff[0] = tr;
  } // Outer loop
  return iters;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int update() {
  int iters = 0;
  Real final_rsq, CKU0 = kappa * clov_c / (u0 * u0 * u0);
  Real cg_time[2], old_cg_time[2], next_cg_time[2];
  double starttrlogA = 0.0;
#ifdef HMC_ALGORITHM
  double startaction = 0, endaction, endtrlogA, change;
  Real xrandom;
#endif
  double gnorm = 0.0, fnorm[2];
  fnorm[0] = 0.0;
  fnorm[1] = 0.0;
  max_gf = 0.0;
  max_ff[0] = 0.0;
  max_ff[1] = 0.0;

#ifdef NHYP_JACOBI
  int i;
  for (i = 0; i < JACOBI_HIST_MAX; i++)
    jacobi_hist[i] = 0;
#endif

#ifndef PHI_ALGORITHM   // Require Phi algorithm
  node0_printf("Only works for phi algorithm\n");
  exit(1);
#endif

  // Refresh the momenta
  ranmom();

  // DIMFxDIMF link created from NCOLxNCOL linkf after each update,
  // then DIMF gauge field is switched to antiperiodic BCs in fermion_rep()
  // Initial DIMF links are set in update()
  // See also sf_make_boundary.c

  // Generate a pseudofermion configuration only at start
  make_clov(CKU0);
#ifdef LU
  starttrlogA = make_clovinv(ODD);
#endif
  grsource_w();
  old_cg_time[0] = -1;
  old_cg_time[1] = -1;
  cg_time[0] = -1;
  cg_time[1] = -1;

  // Do CG to get both psi = (M^dag M)^(-1) chi
  iters += congrad_cl_m(niter, rsqmin, &final_rsq,
                        F_OFFSET(chi[0]), F_OFFSET(psi[0]), 0.0);
  if (num_masses == 2) {
    iters += congrad_cl_m(niter, rsqmin, &final_rsq,
                          F_OFFSET(chi[1]), F_OFFSET(psi[1]), shift);
  }
#ifdef CG_DEBUG
  checkmul();
#endif

  cg_time[0] = 0.0;
  cg_time[1] = 0.0;

#ifdef HMC_ALGORITHM    // Find action
  startaction = d_action();
  startaction -= 2.0 * starttrlogA;
#ifdef CG_DEBUG
  node0_printf("startaction = %g\n", startaction);
#endif

  // Copy link field to old_link
  gauge_field_copy_f(F_OFFSET(linkf[0]), F_OFFSET(old_linkf[0]));
#endif

  // Do microcanonical updating
  iters += update_step(old_cg_time, cg_time, next_cg_time, fnorm, &gnorm);

#ifdef HMC_ALGORITHM    // Find action, then accept or reject
#ifdef LU
  // update_step() provides returntrlogA for us to reuse
  endtrlogA = returntrlogA;
#endif
  // Do CG to get both psi = (M^dag M)^(-1) chi
  iters += congrad_cl_m(niter, rsqmin, &final_rsq,
                         F_OFFSET(chi[0]), F_OFFSET(psi[0]), 0.0);
  if (num_masses == 2) {
    iters += congrad_cl_m(niter, rsqmin, &final_rsq,
                          F_OFFSET(chi[1]), F_OFFSET(psi[1]), shift);
  }

  endaction = d_action();
  endaction -= (double)2.0 * endtrlogA;
#ifdef CG_DEBUG
  node0_printf("endaction = %g\n", endaction);
#endif

  change = endaction - startaction;

  // Reject configurations giving overflow
  if (fabs((double)change)>1e20) {
    node0_printf("WARNING: Correcting apparent overflow: Delta S = %.4g\n",
                 change);
    change = 1.0e20;
  }

  // Decide whether to accept.  If not, copy old link field back
  // Careful: must generate only one random number for whole lattice
  if (this_node == 0)
    xrandom = myrand(&node_prn);

  broadcast_float(&xrandom);
    if (exp(-change) < (double)xrandom) {
      gauge_field_copy_f(F_OFFSET(old_linkf[0]), F_OFFSET(linkf[0]));
      fermion_rep();
      node0_printf("REJECT: ");
    }
    else
      node0_printf("ACCEPT: ");

  node0_printf("delta S = %.4g start %.12g end %.12g\n",
               change, startaction, endaction);
#endif
  free_clov();    // Needed for Phi algorithm in addition to HMC

  node0_printf("MONITOR_FORCE_FERMION0 %.4g %.4g\n",
               fnorm[0]/(double)(2 * nsteps[0]), max_ff[0]);
  node0_printf("MONITOR_FORCE_FERMION1 %.4g %.4g\n",
               fnorm[1]/(double)(4 * nsteps[0] * nsteps[1]), max_ff[1]);
  // gnorm divided by nsteps_gauge every time gauge_update_step called
  node0_printf("MONITOR_FORCE_GAUGE %.4g %.4g\n",
               gnorm/(double)(4 * nsteps[0] * nsteps[1]), max_gf);

#ifdef NHYP_JACOBI
  jacobi_total = 0;
  jacobi_avrg = 0.0;
  for (i = 0; i < JACOBI_HIST_MAX; i++) {
    jacobi_total += jacobi_hist[i];
    jacobi_avrg += (Real)((i + 1) * jacobi_hist[i]);
  }
  jacobi_avrg /= (Real)jacobi_total;
  node0_printf("MONITOR_JACOBI_TOTL %d\n", jacobi_total);
  if (jacobi_total > 0) {
    node0_printf("MONITOR_JACOBI_AVRG %.4g\n", jacobi_avrg);
    node0_printf("MONITOR_JACOBI_HIST   ");
    for (i = 0; i < JACOBI_HIST_MAX; i++)
      node0_printf("%.4g ",  (Real)jacobi_hist[i] / (Real)jacobi_total);
    node0_printf("\n");
  }
#endif
  return iters;
}
// -----------------------------------------------------------------
