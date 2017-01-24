#include "cl_dyn_includes.h"

/* two-mass version. if num_masses = 1, then usual code. If 2, then
chi[0] and psi[0] are the ``slow'' (no-shift) inverter and chi[1] and psi[1]
are the fast (shifted) pseudofermions */

// Construct a gaussian random vector, g_rand, and chi = Mdag.g_rand
// Also clear psi, since zero is our best guess with a new random source

/* Use the LU preconditioned fermion matrix, where
  the fermion spinors live on even sites only.  In other words, if
  Dslash_oe is the dslash operator with its source on even sites and
  its result on odd sites, etc.:
   M = A_even - kappa^2 * Dslash_eo * (A_odd)^{-1} * Dslash_oe
*/

void grsource() {
  register int i, j, k;
  register site *s;
  int index = num_masses - 1, iters = 0;
  wilson_vector twvec;
#ifdef DEBUG_CHECK
  double rr = 0.0;
#endif

  /* First the original field or the simply shifted one in the 2PF case */
  // index is 0 for num_masses = 1 or 1 for num_masses = 2
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
        psi[index][i].d[k].c[j] = cmplx(0.0, 0.0);
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

  // chi[index] <-- Mdag g_rand and the num_masses = 1 case is done
  fermion_op(g_rand, chi[index], MINUS, EVEN);
  if (num_masses == 2) {    // Subtract i * shift * gamma5.g_rand
    FOREVENSITES(i, s) {
      mult_by_gamma(&(g_rand[i]), &twvec, GAMMAFIVE);
      c_scalar_mult_sum_wvec(&twvec, &mishift, &(chi[1][i]));
    }

  // Now the slightly fancier case of chi[0] for two pseudo-fermions
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

    // Now invert with the shift and write the result into psi[0]
    iters += congrad(0, shift, EVEN);

    // mp <-- Mtilde psi[0]
    // chi[0] <-- Mdag mp
    fermion_op(psi[0], mp, PLUS, EVEN);
    FOREVENSITES(i, s) {
      // Add i * shift * gamma5.psi[0] to mp
      mult_by_gamma(&(psi[0][i]), &twvec, GAMMAFIVE);
      c_scalar_mult_sum_wvec(&twvec, &ishift, &(mp[i]));

      // Finally add i * shift * gamma5.mp to chi[0]
      mult_by_gamma(&(mp[i]), &twvec, GAMMAFIVE);
      c_scalar_mult_sum_wvec(&twvec, &ishift, &(chi[0][i]));
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Check congrad by comparing chi[level] vs. Mdag.M psi[level]
// On average these should differ by less than sqrt(rsqmin)
// Only warn about fluctuations at least 5x larger than that tolerance
// Must run CG before calling, of course
void checkmul(int level, Real mshift) {
  register int i, j, k;
  register site *s;
  Real MSq = mshift * mshift, tol = 5.0 * sqrt(rsqmin), tr;
  wilson_vector twvec;

#ifdef DEBUG_CHECK
  printf("Running checkmul with mshift %.4g\n", mshift);
#endif
  // Multiply solution psi[level] by Mdag.M
  fermion_op(psi[level], mp, PLUS, EVEN);
  fermion_op(mp, mp, MINUS, EVEN);
  FOREVENSITES(i, s) {
    // Add mshift^2 psi[level] to complete Mdag.M psi[level]
    if (mshift > 0.0)
      scalar_mult_sum_wvec(&(psi[level][i]), MSq, &(mp[i]));

    // Compare to source in chi[level]
#ifdef DEBUG_CHECK
    printf("Site %d %d %d %d\n", s->x, s->y, s->z, s->t);
#endif
    copy_wvec(&(chi[level][i]), &twvec);
    for (k = 0; k < 4; k++) {
      for (j = 0; j < DIMF; j++) {
        tr = twvec.d[k].c[j].real - mp[i].d[k].c[j].real;
        if (fabs(tr) > tol) {
          printf("real %d %d %d: %.4g - %.4g = %.4g > %.4g\n", i, k, j,
                 twvec.d[k].c[j].real, mp[i].d[k].c[j].real, tr, tol);
        }
        tr = twvec.d[k].c[j].imag - mp[i].d[k].c[j].imag;
        if (fabs(tr) > tol) {
            printf("imag %d %d %d: %.4g - %.4g = %.4g > %.4g\n", i, k, j,
                   twvec.d[k].c[j].imag, mp[i].d[k].c[j].imag, tr, tol);
        }
      }
    }
  }
}
// -----------------------------------------------------------------
