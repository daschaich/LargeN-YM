// -----------------------------------------------------------------
// Solve (Mdag.M) psi = chi for clover fermions in dynamical HMC updates
// This version allows an input mshift for Hasenbusch preconditioning

/* Use the LU preconditioned fermion matrix, where
  the fermion spinors live on even sites only.  In other words, if
  Dslash_oe is the dslash operator with its source on even sites and
  its result on odd sites, etc.:
   M = A_e - kappa^2 * Dslash_eo * (A_o)^{-1} * Dslash_oe
*/

// This version looks at the initial vector every niter passes
// The level tells us which source and initial guess/answer to use:
//   chi[level] and psi[level], respectively
// "r" is the residual vector
// "p" and "mp" are working vectors for the conjugate gradient
// niter = maximum number of iterations
// rsqmin = desired rsq, quit when we reach rsq = rsqmin*source_norm.
#include "cl_dyn_includes.h"

int congrad(int niter, Real rsqmin, Real *final_rsq_ptr,
            int level, Real mshift) {

  register int i;
  register site *s;
  int iteration = 0;
  Real a = 0.0, b, kappaSq = -kappa * kappa, MSq = mshift * mshift;
  double rsqstop = 0.0, rsq = 0.0, oldrsq = 0.0, dtime = -dclock();
  double source_norm = 0.0, pkp = 0.0;      // pkp = p.K.p
  msg_tag *tag[8], *tag2[8];

#ifdef TIMING
  TIC(0)
#endif

  // Allocate fields, copy source and initial guess
  wilson_vector *chi = malloc(sites_on_node * sizeof(*chi));
  wilson_vector *psi = malloc(sites_on_node * sizeof(*psi));
  wilson_vector *tmp = malloc(sites_on_node * sizeof(*tmp));
  wilson_vector *mp  = malloc(sites_on_node * sizeof(*mp));
  wilson_vector *p   = malloc(sites_on_node * sizeof(*p));
  wilson_vector *r   = malloc(sites_on_node * sizeof(*r));
  FORALLSITES(i, s) {
    copy_wvec(&(s->chi[level]), &(chi[i]));
    copy_wvec(&(s->psi[level]), &(psi[i]));
  }

start:
  /* mp <-- Mdag.M on psi
     r, p <-- chi - mp
     rsq = |r|^2
     source_norm = |chi|^2
     */
  rsq = 0.0;
  source_norm = 0.0;
  mult_ldu_field(psi, tmp, EVEN);
  dslash_w_field_special(psi, mp, PLUS, ODD, tag, 0);
  mult_ldu_field(mp, tmp, ODD);
  dslash_w_field_special(tmp, mp, PLUS, EVEN, tag2, 0);
  FOREVENSITES(i, s) {
    scalar_mult_wvec(&(mp[i]), kappaSq, &(mp[i]));
    sum_wvec(&(tmp[i]), &(mp[i]));
  }

  // Above applied M on psi, now apply Mdag
  mult_ldu_field(mp, tmp, EVEN);
  dslash_w_field_special(mp, mp, MINUS, ODD, tag, 1);
  mult_ldu_field(mp, tmp, ODD);
  dslash_w_field_special(tmp, mp, MINUS, EVEN, tag2, 1);
  FOREVENSITES(i, s) {
    scalar_mult_wvec(&(mp[i]), kappaSq, &(mp[i]));
    sum_wvec(&(tmp[i]), &(mp[i]));

    // Now mp holds Mdag.M applied to psi
    // Add mshift^2 psi
    if (mshift > 0.0)
      scalar_mult_sum_wvec(&(psi[i]), MSq, &(mp[i]));

    sub_wilson_vector(&(chi[i]), &(mp[i]), &(r[i]));
    p[i] = r[i];
    magsq_wvec_sum(&(chi[i]), &source_norm);
    magsq_wvec_sum(&(r[i]), &rsq);
  }
  g_doublesum(&source_norm);
  g_doublesum(&rsq);
  iteration++ ;         // Count number of Mdag.M multiplications
  total_iters++;

  rsqstop = rsqmin * source_norm;
#ifdef CG_DEBUG
  double mflops;
  node0_printf("CG source_norm = %.4g ", source_norm);
  node0_printf("--> rsqstop = %.4g\n", rsqstop);
  node0_printf("CG iter %d, rsq %.4g, pkp %.4g, a %.4g\n",
               iteration, rsq, pkp, a);
#endif
  if (rsq <= rsqstop) {
    *final_rsq_ptr = (Real)rsq;
    FORALLUPDIR(i) {
      cleanup_gather(tag[i]);
      cleanup_gather(tag2[i]);
      cleanup_gather(tag[OPP_DIR(i)]);
      cleanup_gather(tag2[OPP_DIR(i)]);
    }

    cleanup_dslash_temps();
    cleanup_tmp_links();

    FORALLSITES(i, s)
      copy_wvec(&(psi[i]), &(s->psi[level]));

    free(psi);
    free(chi);
    free(tmp);
    free(mp);
    free(p);
    free(r);

#ifdef TIMING
    TOC(0, time_dcongrad)
#endif

    return iteration;
  }

  // Main loop -- do until convergence or time to restart
  /*
     oldrsq <-- rsq
     mp <-- Mdag.M p
     pkp <-- p Mdag.M p
     a <-- rsq / pkp
     psi <-- psi + a * p
     r <-- r - a * mp
     rsq <-- |r|^2
     b <-- rsq / oldrsq
     p <-- r + b * p
     */
  do {
    oldrsq = rsq;
    pkp = 0.0;
    mult_ldu_field(p, tmp, EVEN);
    dslash_w_field_special(p, mp, PLUS, ODD, tag, 1);
    mult_ldu_field(mp, tmp, ODD);
    dslash_w_field_special(tmp, mp, PLUS, EVEN, tag2, 1);
    FOREVENSITES(i, s) {
      scalar_mult_wvec(&(mp[i]), kappaSq, &(mp[i]));
      sum_wvec(&(tmp[i]), &(mp[i]));
    }

    mult_ldu_field(mp, tmp, EVEN);
    dslash_w_field_special(mp, mp, MINUS, ODD, tag, 1);
    mult_ldu_field(mp, tmp, ODD);
    dslash_w_field_special(tmp, mp, MINUS, EVEN, tag2, 1);
    FOREVENSITES(i, s) {
      scalar_mult_wvec(&(mp[i]), kappaSq, &(mp[i]));
      sum_wvec(&(tmp[i]), &(mp[i]));

      // Add mshift^2 psi
      if (mshift > 0.0)
        scalar_mult_sum_wvec(&(p[i]), MSq, &(mp[i]));
      wvec_rdot_sum(&(p[i]), &(mp[i]), &pkp);
    }
    g_doublesum(&pkp);
    iteration++;
    total_iters++;

    a = (Real)(rsq / pkp);
    rsq = 0.0;
    FOREVENSITES(i, s) {
      scalar_mult_sum_wvec(&(p[i]), a, &(psi[i]));
      scalar_mult_dif_wvec(&(mp[i]), a, &(r[i]));
      magsq_wvec_sum(&(r[i]), &rsq);
    }
    g_doublesum(&rsq);
#ifdef CG_DEBUG
    node0_printf("CG iter %d, rsq %.4g, pkp %.4g, a %.4g\n",
                 iteration, rsq, pkp, a);
#endif

    if (rsq <= rsqstop) {
      *final_rsq_ptr = (Real)rsq;
      FORALLUPDIR(i) {
        cleanup_gather(tag[i]);
        cleanup_gather(tag2[i]);
        cleanup_gather(tag[OPP_DIR(i)]);
        cleanup_gather(tag2[OPP_DIR(i)]);
      }
      dtime += dclock();
#ifdef CG_DEBUG
      mflops = (double)(2840.0 * volume * iteration
                        / (1e6 * dtime * numnodes()));
      node0_printf("CG time = %.4g iters = %d mflops = %.4g\n",
                   dtime, iteration, mflops);
#endif
      cleanup_dslash_temps();
      cleanup_tmp_links();

      FORALLSITES(i, s)
        copy_wvec(&(psi[i]), &(s->psi[level]));
      free(psi);
      free(chi);
      free(tmp);
      free(mp);
      free(p);
      free(r);

#ifdef TIMING
      TOC(0, time_dcongrad)
#endif

      return iteration;
    }

    b = rsq / oldrsq;
    FOREVENSITES(i, s) {
      scalar_mult_wvec(&(p[i]), b, &(p[i]));
      sum_wvec(&(r[i]), &(p[i]));
    }
  } while (iteration % niter != 0);

  FORALLUPDIR(i) {
    cleanup_gather(tag[i]);
    cleanup_gather(tag2[i]);
    cleanup_gather(tag[OPP_DIR(i)]);
    cleanup_gather(tag2[OPP_DIR(i)]);
  }

  // Hard-coded number of restarts in defines.h
  if (iteration < CONGRAD_RESTART * niter)
    goto start;
  *final_rsq_ptr = rsq;
  if (rsq > rsqstop) {
    node0_printf("WARNING: CG did not converge, size_r = %.2g\n",
                 sqrt(rsq / source_norm));
  }

  cleanup_dslash_temps();
  cleanup_tmp_links();

  FORALLSITES(i, s)
    copy_wvec(&(psi[i]), &(s->psi[level]));

  free(psi);
  free(chi);
  free(tmp);
  free(mp);
  free(p);
  free(r);

#ifdef TIMING
  TOC(0, time_dcongrad)
#endif

  return iteration;
}
// -----------------------------------------------------------------
