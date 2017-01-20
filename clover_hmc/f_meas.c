// -----------------------------------------------------------------
/* Measure fermionic observables:
    psi-bar-psi, fermion action, energy and pressure, psi-bar-gamma_5-psi

   In this version, M is NOT the LU-preconditioned matrix
   The preconditioned matrix is MM
   M = A - kappa*(Dslash_eo + Dslash_oe)
   MM = A_e - kappa^2 * Dslash_eo * (A_o)^{-1} * Dslash_oe

   Output is an FMES line containing the expectation value of
   the real part of Tr(1/M),
   the imaginary part of  Tr(1/M),
   the trace of D_slash(t),
   the average (over i)  trace of D_slash(i),
   the fermion action,
   the real part of Tr(gamma_5/M),
   the imaginary part of Tr(gamma_5/M).

   D_slash = (1+gamma_mu)U*delta_+ + (1-gamma_mu)U_adjoint*delta_-

   These must be post-processed into the thermodynamic observables:
   $2 = Re Tr(1/M)
   $3 = Im Tr(1/M)
   $4 = Tr(Dslash_time)
   $5 = Tr(sum_i Dslash_i)/3
   $6 = Re Tr(gamma_5/M)
   $7 = Im Tr(gamma_5/M)
   4 = number of Dirac components
   3 = number of colors

   pbp = (nflavors/2) * 4 * kappa * $2
   entropy = (nflavors/2) * 2 * kappa * (-$4 + $5)
   energy = (nflavors/2) * 2 * kappa * (-$4 + (4*3 - $2)*tderiv)
   pressure = (nflavors/2) * 2 * kappa * ($5 - (4*3 - $2)*3*sderiv)
    where tderiv = partial 1/kappa_c / partial alpha_t
    where sderiv = partial 1/kappa_c / partial alpha_s

   These are the entropy, etc. summed over color and flavor.
*/
#include "cl_dyn_includes.h"

int f_meas() {
  register int i, j, k, dir, iters;
  register site *s;
  double faction = 0.0, dslash_time = 0.0, dslash_space = 0.0;
  double_complex pbp = cmplx(0.0, 0.0), pbg5p = cmplx(0.0, 0.0);
  half_wilson_vector hwv0, hwv1;
  wilson_vector wv0;
  msg_tag *tag0, *tag1;

  // Gaussian random vector
  FORALLSITES(i, s) {
    for (k = 0; k < 4; k++) {
      for (j = 0; j < DIMF; j++) {
#ifdef SITERAND
        s->g_rand.d[k].c[j].real = gaussian_rand_no(&(s->site_prn));
        s->g_rand.d[k].c[j].imag = gaussian_rand_no(&(s->site_prn));
#else
        s->g_rand.d[k].c[j].real = gaussian_rand_no(&node_prn);
        s->g_rand.d[k].c[j].imag = gaussian_rand_no(&node_prn);
#endif
      }
    }
  }

  // Copy gaussian source to chi[0]
  FORALLSITES(i, s)
    copy_wvec(&(s->g_rand), &(s->chi[0]));

  // Load inversion control structure
  qic.start_flag = 0;   // Use zero initial guess for psi

  // Load Dirac matrix parameters, including temporaries
  qic.wv1 = F_OFFSET(tmp);
  qic.wv2 = F_OFFSET(mp);

  // Invert M, put result in psi[0]
  iters = wilson_invert(F_OFFSET(chi[0]), F_OFFSET(psi[0]), F_OFFSET(r),
                        cgilu_cl, &qic, (void *)&dcp);

#ifdef DEBUG_CHECK
  // Multiply by M and see if g_rand is restored
  dslash_w_site(psi[0], mp, PLUS, EVENANDODD);
  FORALLSITES(i, s) {
    scalar_mult_wvec(&(mp[i]), -kappa, &(mp[i]));
    sum_wvec(&(psi[0][i]), &(s->mp[i]));
  }
  FORALLSITES(i, s) {
    for (j = 0; j < 4; j++) {
      for (k = 0; k < 3; k++) {
        tr = g_rand[i].d[j].c[k].real - mp[i].d[j].c[k].real;
        if (fabs(tr) > IMAG_TOL) {
          printf("real %d %d %d: %.4g - %.4g = %.4g > %.4g\n", i, j, k,
                 g_rand[i].d[k].c[j].real, mp[i].d[k].c[j].real,
                 tr, IMAG_TOL);
        }
        tr = g_rand[i].d[j].c[k].imag - mp[i].d[j].c[k].imag;
        if (fabs(tr) > IMAG_TOL) {
          printf("imag %d %d %d: %.4g - %.4g = %.4g > %.4g\n", i, j, k,
                 g_rand[i].d[j].c[k].imag, mp[i].d[j].c[k].imag,
                 tr, IMAG_TOL);
        }
      }
    }
  }
#endif

  // pbp = g_rand.psi
  // pbg5p = g_rand.gamma5.psi
  FORALLSITES(i, s) {
    wvec_dot_sum(&(s->g_rand), &(s->psi[0]), &pbp);
    mult_by_gamma(&(s->psi[0]), &wv0, GAMMAFIVE);
    wvec_dot_sum(&(s->g_rand), &wv0, &pbg5p);
  }

  /* fermion energy and pressure */
  FORALLUPDIR(dir) {
    /* multiply g_rand by one component of Dslash_adjoint, result in p */
    FORALLSITES(i, s) {
      wp_shrink(&(s->g_rand), &(s->htmp[0]), dir, MINUS);
      wp_shrink(&(s->g_rand), &hwv1, dir, PLUS);
      mult_adj_su3_mat_hwvec(&(s->link[dir]), &hwv1, &(s->htmp[1]));
    }
    tag0 = start_gather_site(F_OFFSET(htmp[0]), sizeof(half_wilson_vector),
                             dir, EVENANDODD, gen_pt[0]);
    tag1 = start_gather_site(F_OFFSET(htmp[1]), sizeof(half_wilson_vector),
                             OPP_DIR(dir), EVENANDODD, gen_pt[1]);
    wait_gather(tag0);
    wait_gather(tag1);
    FORALLSITES(i, s) {
      mult_su3_mat_hwvec(&(s->link[dir]), (half_wilson_vector *)(gen_pt[0][i]),
                         &hwv0);
      wp_grow(&hwv0, &(s->p), dir, MINUS);
      wp_grow((half_wilson_vector *)(gen_pt[1][i]), &wv0, dir, PLUS);
      sum_wvec(&wv0, &(s->p));
    }
    cleanup_gather(tag0);
    cleanup_gather(tag1);

    /* dot product with psi, result into energy or pressure */
    FORALLSITES(i, s) {
      if (dir == TUP)
        wvec_rdot_sum(&(s->psi[0]), &(s->p), &dslash_time);
      else
        wvec_rdot_sum(&(s->psi[0]), &(s->p), &dslash_space);
    }
  }
  g_doublesum(&dslash_time);
  g_doublesum(&dslash_space);
  g_complexsum(&pbp);
  g_complexsum(&pbg5p);

  CDIVREAL(pbp, volume, pbp);
  CDIVREAL(pbg5p, volume, pbg5p);
  dslash_time /= (double)volume;
  dslash_space /= (double)(3.0 * volume);
  faction = pbp.real - kappa * (dslash_time + 3.0 * dslash_space);

  // Check that pbg5p is purely real up to machine precision
  if (fabs(pbg5p.imag) > IMAG_TOL)
    node0_printf("WARNING: Im(pbg5p) = %.4g > %.4g\n", pbg5p.imag, IMAG_TOL);

  // Print results
  node0_printf("FMES %.8g %.8g %.8g %.8g %.8g %.8g\n", pbp.real, pbp.imag,
               dslash_time, dslash_space, faction, pbg5p.real);
  fflush(stdout);
  return iters;
}
// -----------------------------------------------------------------
