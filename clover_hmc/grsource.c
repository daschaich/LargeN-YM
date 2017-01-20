#include "cl_dyn_includes.h"

/* two-mass version. if num_masses = 1, then usual code. If 2, then
chi[0] and psi[0] are the ``slow'' (no-shift) inverter and chi[1] and psi[1]
are the fast (shifted) pseudofermions */

// Construct a gaussian random vector, g_rand, and chi = Mdag.g_rand
/* also clear psi, since zero is our best guess for the solution with a
   new random chi field. */

/* Use the LU preconditioned fermion matrix, where
  the fermion spinors live on even sites only.  In other words, if
  Dslash_oe is the dslash operator with its source on even sites and
  its result on odd sites, etc.:
   M = A_even - kappa^2 * Dslash_eo * (A_odd)^{-1} * Dslash_oe
*/

void grsource_w() {
  register int i, j, k;
  register site *s;
  int kind1, iters = 0;
  Real final_rsq;
  double kappaSq = -kappa * kappa;
  complex ishift = cmplx(0.0, shift), mishift = cmplx(0.0, -shift);
  wilson_vector twvec, twvec2;
#ifdef DEBUG_CHECK
  double rr = 0.0;
#endif

  /* First the original field or the simply shifted one in the 2PF case */
  if (num_masses == 1)
    kind1 = 0;
  else
    kind1 = 1;

  FOREVENSITES(i, s) {
    for (k = 0; k < 4; k++) {
      for (j = 0; j < DIMF; j++) {
#ifdef SITERAND
        s->g_rand.d[k].c[j].real = gaussian_rand_no(&(s->site_prn));
        s->g_rand.d[k].c[j].imag = gaussian_rand_no(&(s->site_prn));
#else
        s->g_rand.d[k].c[j].real = gaussian_rand_no(&node_prn);
        s->g_rand.d[k].c[j].imag = gaussian_rand_no(&node_prn);
#endif
        s->psi[kind1].d[k].c[j] = cmplx(0.0, 0.0);
      }
    }
#ifdef DEBUG_CHECK
    wvec_rdot_sum(&(s->g_rand), &(s->g_rand), &rr);
#endif
  }
#ifdef DEBUG_CHECK
  g_doublesum(&rr);
  node0_printf("rr %.4g\n", rr);
#endif

  // chi <-- Mdag g_rand
  mult_ldu_site(F_OFFSET(g_rand), F_OFFSET(tmp), EVEN);
  dslash_w_site(F_OFFSET(g_rand), F_OFFSET(tmp), MINUS, ODD);
  mult_ldu_site(F_OFFSET(tmp), F_OFFSET(g_rand), ODD);
  dslash_w_site(F_OFFSET(g_rand), F_OFFSET(chi[kind1]), MINUS, EVEN);
  FOREVENSITES(i, s) {
    if (num_masses == 2) {
      scalar_mult_add_wvec(&(s->tmp), &(s->chi[kind1]), kappaSq, &twvec);
      /* That was Mdag, now we need to subtract -i*shift*gamma5*p */
      mult_by_gamma(&(s->g_rand), &twvec2, GAMMAFIVE);
      c_scalar_mult_add_wvec(&twvec, &twvec2, &mishift, &(s->chi[kind1]));
    }
    else {
      scalar_mult_wvec(&(s->chi[kind1]), kappaSq, &(s->chi[kind1]));
      sum_wvec(&(s->tmp), &(s->chi[kind1]));
    }
  }

  /* Now the slightly more fancy case of the 2nd PF field: chi[0] and psi[0] */
  if (num_masses == 2) {
#ifdef DEBUG_CHECK
    node0_printf("CALL 2nd PF\n");
    rr = 0.0;
#endif
    /* The first contribution is just the random number itself
       so we might as well generate it in chi[0] directly */
    FOREVENSITES(i, s) {
      for (k = 0; k < 4; k++) {
        for (j = 0; j < DIMF; j++) {
          /* Comment these lines out, set s->chi[0].d[k].c[j] = cmplx(1.0, 0.0);
             and set ``step = 0.0'' in the input file and the 1 and 2 PF codes
             should produce identical results */
#ifdef SITERAND
          s->chi[0].d[k].c[j].real = gaussian_rand_no(&(s->site_prn));
          s->chi[0].d[k].c[j].imag = gaussian_rand_no(&(s->site_prn));
#else
          s->chi[0].d[k].c[j].real = gaussian_rand_no(&node_prn);
          s->chi[0].d[k].c[j].imag = gaussian_rand_no(&node_prn);
#endif
          s->psi[0].d[k].c[j] = cmplx(0.0, 0.0);
        }
      }
#ifdef DEBUG_CHECK
      wvec_rdot_sum(&(s->chi[0]), &(s->chi[0]), &rr);
#endif
    }
#ifdef DEBUG_CHECK
    g_doublesum(&rr);
    node0_printf("rr %.4g\n", rr);
#endif

    /* But now something fancy has to happen. First we invert with the shift
       and write the result into psi[0] */
    iters += congrad(niter, rsqmin, &final_rsq,
                     F_OFFSET(chi[0]), F_OFFSET(psi[0]), shift);

    // chi <-- Mtilde psi1
    // chi <-- Mdag chi[0]
    mult_ldu_site(F_OFFSET(psi[0]), F_OFFSET(tmp), EVEN);
    dslash_w_site(F_OFFSET(psi[0]), F_OFFSET(tmp), PLUS, ODD);
    mult_ldu_site(F_OFFSET(tmp), F_OFFSET(psi[0]), ODD);
    dslash_w_site(F_OFFSET(psi[0]), F_OFFSET(p), PLUS, EVEN);
    FOREVENSITES(i, s) {
      scalar_mult_add_wvec(&(s->tmp), &(s->p), kappaSq, &twvec);

      /* That was M, now we need to add i*shift*gamma5*p */
      mult_by_gamma(&(s->psi[0]), &twvec2, GAMMAFIVE);
      c_scalar_mult_sum_wvec(&twvec2, &ishift, &twvec);
      /* And the final step is to multiply by i*shift*gamma5
         before adding it to chi[0] */
      mult_by_gamma(&twvec, &twvec2, GAMMAFIVE);
      c_scalar_mult_sum_wvec(&twvec2, &ishift, &(s->chi[0]));
    }
  }
}


/* Check congrad by multiplying sol by Mdag.M, compare result to src */
// Must run CG before calling
void checkmul(field_offset src, field_offset sol, Real mshift) {
  register int i, j, k;
  register site *s;
  Real kappaSq = -kappa * kappa, MSq = mshift * mshift, tr;
  wilson_vector wvec;

#ifdef DEBUG_CHECK
  printf("CHECKMUL: starting with mshift %.4g\n", mshift);
#endif
  // Multiply by Mdag.M
  mult_ldu_site(sol, F_OFFSET(tmp), EVEN);
  dslash_w_site(sol, sol, PLUS, ODD);
  mult_ldu_site(sol, F_OFFSET(mp), ODD);
  dslash_w_site(F_OFFSET(mp), F_OFFSET(mp), PLUS, EVEN);
  FOREVENSITES(i, s) {
    scalar_mult_wvec(&(s->mp), kappaSq, &(s->mp));
    sum_wvec(&(s->tmp), &(s->mp));
  }
  mult_ldu_site(F_OFFSET(mp), F_OFFSET(tmp), EVEN);
  dslash_w_site(F_OFFSET(mp), F_OFFSET(mp), MINUS, ODD);
  mult_ldu_site(F_OFFSET(mp), F_OFFSET(tmp), ODD);
  dslash_w_site(F_OFFSET(tmp), F_OFFSET(p), MINUS, EVEN);
  FOREVENSITES(i, s) {
    scalar_mult_wvec(&(s->p), kappaSq, &(s->p));
    sum_wvec(&(s->tmp), &(s->p));
    scalar_mult_sum_wvec((wilson_vector *)F_PT(s, sol), MSq, &(s->p));
  }

  // Compare to source
  FOREVENSITES(i, s) {
#ifdef DEBUG_CHECK
    printf("Site %d %d %d %d\n", s->x, s->y, s->z, s->t);
#endif
    copy_wvec((wilson_vector *)F_PT(s, src), &wvec);
    for (k = 0; k < 4; k++) {
      for (j = 0; j < DIMF; j++) {
        tr = wvec.d[k].c[j].real - s->p.d[k].c[j].real;
        if (fabs(tr) > IMAG_TOL) {
          printf("real %d %d %d: %.4g - %.4g = %.4g > %.4g\n", i, k, j,
                 wvec.d[k].c[j].real, s->p.d[k].c[j].real, tr, IMAG_TOL);
        }
        tr = wvec.d[k].c[j].imag - s->p.d[k].c[j].imag;
        if (fabs(tr) > IMAG_TOL) {
            printf("imag %d %d %d: %.4g - %.4g = %.4g > %.4g\n", i, k, j,
                   wvec.d[k].c[j].imag, s->p.d[k].c[j].imag, tr, IMAG_TOL);
        }
      }
    }
  }
}
// -----------------------------------------------------------------
