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
  Real mkappaSq = -kappa * kappa;
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
        g_rand[i].d[k].c[j].real = gaussian_rand_no(&(s->site_prn));
        g_rand[i].d[k].c[j].imag = gaussian_rand_no(&(s->site_prn));
#else
        g_rand[i].d[k].c[j].real = gaussian_rand_no(&node_prn);
        g_rand[i].d[k].c[j].imag = gaussian_rand_no(&node_prn);
#endif
        psi[kind1][i].d[k].c[j] = cmplx(0.0, 0.0);
      }
    }
#ifdef DEBUG_CHECK
    wvec_rdot_sum(&(g_rand[i]), &(g_rand[i]), &rr);
#endif
  }
#ifdef DEBUG_CHECK
  g_doublesum(&rr);
  node0_printf("rr %.4g\n", rr);
#endif

  // chi <-- Mdag g_rand
  mult_ldu_field(g_rand, tempwvec, EVEN);
  dslash_w_field(g_rand, tempwvec, MINUS, ODD);
  mult_ldu_field(tempwvec, g_rand, ODD);
  dslash_w_field(g_rand, chi[kind1], MINUS, EVEN);
  FOREVENSITES(i, s) {
    if (num_masses == 2) {
      scalar_mult_add_wvec(&(tempwvec[i]), &(chi[kind1][i]), mkappaSq, &twvec);

      /* That was Mdag, now we need to subtract -i*shift*gamma5*p */
      mult_by_gamma(&(g_rand[i]), &twvec2, GAMMAFIVE);
      c_scalar_mult_add_wvec(&twvec, &twvec2, &mishift, &(chi[kind1][i]));
    }
    else {
      scalar_mult_wvec(&(chi[kind1][i]), mkappaSq, &(chi[kind1][i]));
      sum_wvec(&(tempwvec[i]), &(chi[kind1][i]));
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
          /* Comment these lines out, set chi[0][i].d[k].c[j] = cmplx(1.0, 0.0);
             and set ``step = 0.0'' in the input file and the 1 and 2 PF codes
             should produce identical results */
#ifdef SITERAND
          chi[0][i].d[k].c[j].real = gaussian_rand_no(&(s->site_prn));
          chi[0][i].d[k].c[j].imag = gaussian_rand_no(&(s->site_prn));
#else
          chi[0][i].d[k].c[j].real = gaussian_rand_no(&node_prn);
          chi[0][i].d[k].c[j].imag = gaussian_rand_no(&node_prn);
#endif
          psi[0][i].d[k].c[j] = cmplx(0.0, 0.0);
        }
      }
#ifdef DEBUG_CHECK
      wvec_rdot_sum(&(chi[0][i]), &(chi[0][i]), &rr);
#endif
    }
#ifdef DEBUG_CHECK
    g_doublesum(&rr);
    node0_printf("rr %.4g\n", rr);
#endif

    /* But now something fancy has to happen. First we invert with the shift
       and write the result into psi[0] */
    iters += congrad(niter, rsqmin, 0, shift);

    // chi <-- Mtilde psi1
    // chi <-- Mdag chi[0]
    mult_ldu_field(psi[0], tempwvec, EVEN);
    dslash_w_field(psi[0], tempwvec, PLUS, ODD);
    mult_ldu_field(tempwvec, psi[0], ODD);
    dslash_w_field(psi[0], p, PLUS, EVEN);
    FOREVENSITES(i, s) {
      scalar_mult_add_wvec(&(tempwvec[i]), &(p[i]), mkappaSq, &twvec);

      /* That was M, now we need to add i*shift*gamma5*p */
      mult_by_gamma(&(psi[0][i]), &twvec2, GAMMAFIVE);
      c_scalar_mult_sum_wvec(&twvec2, &ishift, &twvec);
      /* And the final step is to multiply by i*shift*gamma5
         before adding it to chi[0] */
      mult_by_gamma(&twvec, &twvec2, GAMMAFIVE);
      c_scalar_mult_sum_wvec(&twvec2, &ishift, &(chi[0][i]));
    }
  }
}


/* Check congrad by multiplying sol by Mdag.M, compare result to src */
// Must run CG before calling
void checkmul(wilson_vector *src, wilson_vector *sol, Real mshift) {
  register int i, j, k;
  register site *s;
  Real mkappaSq = -kappa * kappa, MSq = mshift * mshift, tr;
  wilson_vector wvec;

#ifdef DEBUG_CHECK
  printf("CHECKMUL: starting with mshift %.4g\n", mshift);
#endif
  // Multiply by Mdag.M
  mult_ldu_field(sol, tempwvec, EVEN);
  dslash_w_field(sol, sol, PLUS, ODD);
  mult_ldu_field(sol, mp, ODD);
  dslash_w_field(mp, mp, PLUS, EVEN);
  FOREVENSITES(i, s) {
    scalar_mult_wvec(&(mp[i]), mkappaSq, &(mp[i]));
    sum_wvec(&(tempwvec[i]), &(mp[i]));
  }
  mult_ldu_field(mp, tempwvec, EVEN);
  dslash_w_field(mp, mp, MINUS, ODD);
  mult_ldu_field(mp, tempwvec, ODD);
  dslash_w_field(tempwvec, p, MINUS, EVEN);
  FOREVENSITES(i, s) {
    scalar_mult_wvec(&(p[i]), mkappaSq, &(p[i]));
    sum_wvec(&(tempwvec[i]), &(p[i]));
    scalar_mult_sum_wvec(&(sol[i]), MSq, &(p[i]));
  }

  // Compare to source
  FOREVENSITES(i, s) {
#ifdef DEBUG_CHECK
    printf("Site %d %d %d %d\n", s->x, s->y, s->z, s->t);
#endif
    copy_wvec(&(src[i]), &wvec);
    for (k = 0; k < 4; k++) {
      for (j = 0; j < DIMF; j++) {
        tr = wvec.d[k].c[j].real - p[i].d[k].c[j].real;
        if (fabs(tr) > IMAG_TOL) {
          printf("real %d %d %d: %.4g - %.4g = %.4g > %.4g\n", i, k, j,
                 wvec.d[k].c[j].real, p[i].d[k].c[j].real, tr, IMAG_TOL);
        }
        tr = wvec.d[k].c[j].imag - p[i].d[k].c[j].imag;
        if (fabs(tr) > IMAG_TOL) {
            printf("imag %d %d %d: %.4g - %.4g = %.4g > %.4g\n", i, k, j,
                   wvec.d[k].c[j].imag, p[i].d[k].c[j].imag, tr, IMAG_TOL);
        }
      }
    }
  }
}
// -----------------------------------------------------------------
