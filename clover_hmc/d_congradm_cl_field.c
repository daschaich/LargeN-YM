// -----------------------------------------------------------------
// Solve (Mdag.M) psi = chi for clover fermions in dynamical HMC updates

/* Use the LU preconditioned fermion matrix, where
  the fermion spinors live on even sites only.  In other words, if
  Dslash_oe is the dslash operator with its source on even sites and
  its result on odd sites, etc.:
   M = A_e - kappa^2 * Dslash_eo * (A_o)^{-1} * Dslash_oe
*/

/* adapted to "field" instead of "site" */

// This version looks at the initial vector every "niter" passes
// "chi" is the source vector
// "psi" is the initial guess and answer
// "r" is the residual vector
// "p" and "mp" are working vectors for the conjugate gradient
// niter = maximum number of iterations
// rsqmin = desired rsq, quit when we reach rsq = rsqmin*source_norm.
#include "cl_dyn_includes.h"

int congrad_cl_m(int niter, Real rsqmin, Real *final_rsq_ptr,
                 field_offset src, field_offset dest, Real mshift) {

  register int i;
  register site *s;
  int iteration = 0;
  Real a = 0.0, b, kappaSq = -kappa * kappa, MSq = mshift * mshift;
  double rsqstop = 0.0, rsq = 0.0, oldrsq = 0.0, dtime = -dclock();
  double source_norm = 0.0, pkp = 0.0;      // pkp = p.K.p
  msg_tag *tag[8], *tag2[8];
  void dslash_w_field_special();

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
    copy_wvec((wilson_vector *)F_PT(s, src), &(chi[i]));
    copy_wvec((wilson_vector *)F_PT(s, dest), &(psi[i]));
  }

start:
  /* mp <--  Mdag.M on psi
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
    scalar_mult_wvec(mp + i, kappaSq, mp + i);
    sum_wvec(tmp + i, mp + i);
  }

  // Above applied M on psi, now apply Mdag
  mult_ldu_field(mp, tmp, EVEN);
  dslash_w_field_special(mp, mp, MINUS, ODD, tag, 1);
  mult_ldu_field(mp, tmp, ODD);
  dslash_w_field_special(tmp, mp, MINUS, EVEN, tag2, 1);
  FOREVENSITES(i, s) {
    scalar_mult_wvec(mp + i, kappaSq, mp + i);
    sum_wvec(tmp + i, mp + i);

    // Now mp holds Mdag.M applied to psi
    // Add mshift^2 psi
    if (mshift > 0.0)
      scalar_mult_sum_wvec(psi + i, MSq, mp + i);

    sub_wilson_vector(chi + i, mp + i, r + i);
    p[i] = r[i];
    source_norm += (double)magsq_wvec(chi + i);
    rsq += (double)magsq_wvec(r + i);
  }
  g_doublesum(&source_norm);
  g_doublesum(&rsq);
  iteration++ ;         // Count number of Mdag.M multiplications
  total_iters++;
#ifdef CG_DEBUG
  node0_printf("CG source_norm = %.4g\n", source_norm);
  node0_printf("CG iter %d, rsq %.4g, pkp %.4g, a %.4g\n",
               iteration, rsq, pkp, a);
#endif

  rsqstop = rsqmin * source_norm;
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
      copy_wvec(&(psi[i]), (wilson_vector *)F_PT(s, dest));

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
     oldrsq <- rsq
     mp <- M_adjoint*M*p
     pkp <- p.M_adjoint*M.p
     a <- rsq/pkp
     psi <- psi + a*p
     r <- r - a*mp
     rsq <- |r|^2
     b <- rsq/oldrsq
     p <- r + b*p
     */
  do {
    oldrsq = rsq;
    pkp = 0.0;
    mult_ldu_field(p, tmp, EVEN);
    dslash_w_field_special(p, mp, PLUS, ODD, tag, 1);
    mult_ldu_field(mp, tmp, ODD);
    dslash_w_field_special(tmp, mp, PLUS, EVEN, tag2, 1);
    FOREVENSITES(i, s)
      scalar_mult_add_wvec(tmp + i, mp + i, kappaSq, mp + i);

    mult_ldu_field(mp, tmp, EVEN);
    dslash_w_field_special(mp, mp, MINUS, ODD, tag, 1);
    mult_ldu_field(mp, tmp, ODD);
    dslash_w_field_special(tmp, mp, MINUS, EVEN, tag2, 1);
    FOREVENSITES(i, s) {
      scalar_mult_add_wvec(tmp + i, mp + i, kappaSq, mp + i);

      // Add mshift^2 psi
      if (mshift > 0.0)
        scalar_mult_sum_wvec(p + i, MSq, mp + i);
      pkp += (double)wvec_rdot(p + i, mp + i);
    }
    g_doublesum(&pkp);
    iteration++;
    total_iters++;

    a = (Real)(rsq / pkp);
    rsq = 0.0;
    FOREVENSITES(i, s) {
      scalar_mult_add_wvec(psi + i, p + i, a, psi + i);
      scalar_mult_sum_wvec(mp + i, -a, r + i);
      rsq += (double)magsq_wvec(r + i);
    }
    g_doublesum(&rsq);
    /* node0_printf("congrad2: iter %d, rsq %e, pkp %e, a %e\n",
       iteration,(double)rsq,(double)pkp,(double)a);  */

    if (rsq <= rsqstop) {
      *final_rsq_ptr = (Real)rsq;
      FORALLUPDIR(i) {
        cleanup_gather(tag[i]);
        cleanup_gather(tag2[i]);
        cleanup_gather(tag[OPP_DIR(i)]);
        cleanup_gather(tag2[OPP_DIR(i)]);
      }
      dtime += dclock();
      /* node0_printf("CONGRAD2: time = %e iters = %d mflops = %e\n",
         dtime, iteration,(double)(2840.0*volume*iteration/(1.0e6*dtime*numnodes()))); */
      cleanup_dslash_temps();
      cleanup_tmp_links();

      FORALLSITES(i, s)
        copy_wvec(&(psi[i]), (wilson_vector *)F_PT(s, dest));
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

    b = (Real)(rsq/oldrsq);
    FOREVENSITES(i, s)
      scalar_mult_add_wvec(r + i, p + i, b, p + i);
  } while (iteration % niter != 0);

  FORALLUPDIR(i) {
    cleanup_gather(tag[i]);
    cleanup_gather(tag2[i]);
    cleanup_gather(tag[OPP_DIR(i)]);
    cleanup_gather(tag2[OPP_DIR(i)]);
  }

  /* hard-coded number of restarts in defines.h...  */
  if (iteration < CONGRAD_RESTART * niter)
    goto start;
  *final_rsq_ptr = (Real)rsq;
  if (rsq > rsqstop) {
    node0_printf("WARNING: CG did not converge, size_r = %.2g\n",
                 sqrt(rsq / source_norm));
  }

  cleanup_dslash_temps();
  cleanup_tmp_links();

  FORALLSITES(i, s)
    copy_wvec(&(psi[i]), (wilson_vector *)F_PT(s, dest));

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
