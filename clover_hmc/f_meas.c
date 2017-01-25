// -----------------------------------------------------------------
// Measure fermionic observables: psibar.psi, psibar.gamma5.psi,
//   fermion action, energy and pressure

// Use the LU preconditioned fermion matrix
// Solution and result live on even sites only
// If dslash_oe takes a source on even sites to a result on odd sites, then
//   M = R_e - kappa^2 * dslash_eo * R_o^(-1) * dslash_oe

// Schematically,
//   dslash_oe = (1 + isign * gamma_mu).U_mu(x) src_e(x + mu)
//             + (1 - isign *gamma_mu).Udag_mu(x - mu) src_e(x - mu)
// with a sum over mu and site x in the odd sublattice

// Output is an FMES line containing the following volume averages:
// 1) the real part of Tr[M^(-1)]
// 2) the imag part of Tr[M^(-1)]
// 3) the trace of dslash_t
// 4) the trace of dslash_i averaged over spatial directions i
// 5) the fermion action, which should be 4N on average
// 6) the real part of Tr(gamma5.M^(-1)),
// We check that the imaginary part of Tr(gamma5.M^(-1)) vanishes

// The traces of dslash require solutions on both EVEN and ODD sublattices

// The output above must be post-processed
// into the thermodynamic observables summed over color and flavor:
/* $2 = ReTr[M^(-1)]
   $3 = ImTr[M^(-1)]
   $4 = Tr[dslash_time]
   $5 = Tr[sum_i Dslash_i]/3

   pbp = (nflavors/2) * 4 * kappa * $2
   entropy = (nflavors/2) * 2 * kappa * (-$4 + $5)
   energy = (nflavors/2) * 2 * kappa * (-$4 + (4 * NCOL - $2) * tderiv)
   pressure = (nflavors/2) * 2 * kappa * ($5 - (4 * NCOL - $2) * 3 * sderiv)
    where tderiv = partial 1/kappa_c / partial alpha_t
    where sderiv = partial 1/kappa_c / partial alpha_s
*/
#include "cl_dyn_includes.h"

int f_meas() {
  register int i, j, k, iters;//, dir
  register site *s;
  double faction = 0.0;//, dslash_time = 0.0, dslash_space = 0.0, check = 0.0;
  double_complex pbp = cmplx(0.0, 0.0), pbg5p = cmplx(0.0, 0.0);
//  half_wilson_vector hvec;
  wilson_vector twvec;
//  msg_tag *tag0, *tag1;

  // Gaussian random vector -- simpler than two-mass case in grsource.c
  // Also clear psi[0], since zero is our best guess with a new random source
  FORALLSITES(i, s) {
    for (k = 0; k < 4; k++) {
      for (j = 0; j < DIMF; j++) {
#ifdef SITERAND
        g_rand[i].d[k].c[j].real = gaussian_rand_no(&(s->site_prn));
        g_rand[i].d[k].c[j].imag = gaussian_rand_no(&(s->site_prn));
#else
        g_rand[i].d[k].c[j].real = gaussian_rand_no(&node_prn);
        g_rand[i].d[k].c[j].imag = gaussian_rand_no(&node_prn);
#endif
        psi[0][i].d[k].c[j] = cmplx(0.0, 0.0);
      }
    }
  }

  // Invert on even sites
  // make_clov computes R = 1 - i * CKU0 sigma_{\mu\nu} F_{\mu\nu}
  // make_clovinv inverts inverts R_o on odd sites only
  make_clov(CKU0);
  make_clovinv(ODD);

  // chi[0] <-- Mdag g_rand
  fermion_op(g_rand, chi[0], MINUS, EVEN);

  // Invert Mdag.M with source chi[0], put result in psi[0]
  iters = congrad(0, 0.0, EVEN);

  // Repeat the steps above for odd sites
  free_clov();
  make_clov(CKU0);
  make_clovinv(EVEN);
  fermion_op(g_rand, chi[0], MINUS, ODD);
  iters += congrad(0, 0.0, ODD);
  free_clov();

#ifdef DEBUG_CHECK
  // Hit psi[0] with M and see if g_rand is restored
  // On average these should differ by less than sqrt(rsqmin)
  // Only warn about fluctuations at least 5x larger than that tolerance
  Real tr, tol = 5.0 * sqrt(rsqmin);

  make_clov(CKU0);
  make_clovinv(ODD);
  fermion_op(psi[0], mp, PLUS, EVEN);
  free_clov();
  make_clov(CKU0);
  make_clovinv(EVEN);
  fermion_op(psi[0], mp, PLUS, ODD);
  free_clov();

  FORALLSITES(i, s) {
    for (j = 0; j < 4; j++) {
      for (k = 0; k < NCOL; k++) {
        tr = g_rand[i].d[j].c[k].real - mp[i].d[j].c[k].real;
        if (fabs(tr) > tol) {
          printf("real %d %d %d: %.4g - %.4g = %.4g > %.4g\n", i, j, k,
                 g_rand[i].d[k].c[j].real, mp[i].d[k].c[j].real,
                 tr, tol);
        }
        tr = g_rand[i].d[j].c[k].imag - mp[i].d[j].c[k].imag;
        if (fabs(tr) > tol) {
          printf("imag %d %d %d: %.4g - %.4g = %.4g > %.4g\n", i, j, k,
                 g_rand[i].d[j].c[k].imag, mp[i].d[j].c[k].imag,
                 tr, tol);
        }
      }
    }
  }
#endif

  // pbp = g_rand^dag.psi[0] = g_rand^dag.M^(-1).g_rand
  // pbg5p = g_rand^dag.gamma5.psi[0] = g_rand.gamma5.M^(-1).g_rand
  // faction = chi[0]^dag.psi = g_rand.g_rand
  FORALLSITES(i, s) {
    wvec_dot_sum(&(g_rand[i]), &(psi[0][i]), &pbp);

    mult_by_gamma(&(psi[0][i]), &twvec, GAMMAFIVE);
    wvec_dot_sum(&(g_rand[i]), &twvec, &pbg5p);

    wvec_rdot_sum(&(chi[0][i]), &(psi[0][i]), &faction);
#ifdef DEBUG_CHECK
    wvec_rdot_sum(&(g_rand[i]), &(g_rand[i]), &check);
#endif
  }
  g_dcomplexsum(&pbp);
  g_dcomplexsum(&pbg5p);
  g_doublesum(&faction);
  CMULREAL(pbp, one_ov_vol, pbp);
  CMULREAL(pbg5p, one_ov_vol, pbg5p);
  faction *= one_ov_vol;
#ifdef DEBUG_CHECK
  g_doublesum(&check);
  check *= one_ov_vol;
  node0_printf("CHECK faction: %.4g - %.4g = %.4g\n",
               faction, check, faction - check);
#endif

  // TODO: dslash_time and dslash_space change significantly
  //       when using congrad rather than wilson_invert
  //       However, even in commit 3f59e86 the 'faction' computed below
  //       differed from Tr[gdag.g] by ~10%...
  //       Comment out for now
#if 0
  /* fermion energy and pressure */
  FORALLUPDIR(dir) {
    /* multiply g_rand by one component of Dslash_adjoint, result in p */
    FORALLSITES(i, s) {
      wp_shrink(&(g_rand[i]), &(htmp[0][i]), dir, MINUS);
      wp_shrink(&(g_rand[i]), &hvec, dir, PLUS);
      mult_adj_mat_hwvec(&(s->link[dir]), &hvec, &(htmp[1][i]));
    }
    tag0 = start_gather_field(htmp[0], sizeof(half_wilson_vector),
                              dir, EVENANDODD, gen_pt[0]);
    tag1 = start_gather_field(htmp[1], sizeof(half_wilson_vector),
                              OPP_DIR(dir), EVENANDODD, gen_pt[1]);
    wait_gather(tag0);
    wait_gather(tag1);
    FORALLSITES(i, s) {
      mult_mat_hwvec(&(s->link[dir]), (half_wilson_vector *)(gen_pt[0][i]),
                         &hvec);
      wp_grow(&hvec, &(p[i]), dir, MINUS);
      wp_grow((half_wilson_vector *)(gen_pt[1][i]), &twvec, dir, PLUS);
      sum_wvec(&twvec, &(p[i]));
    }
    cleanup_gather(tag0);
    cleanup_gather(tag1);

    /* dot product with psi[0], result into energy or pressure */
    FORALLSITES(i, s) {
      if (dir == TUP)
        wvec_rdot_sum(&(psi[0][i]), &(p[i]), &dslash_time);
      else
        wvec_rdot_sum(&(psi[0][i]), &(p[i]), &dslash_space);
    }
  }
  g_doublesum(&dslash_time);
  g_doublesum(&dslash_space);
  dslash_time *= one_ov_vol;
  dslash_space *= one_ov_vol;
  check = pbp.real - kappa * (dslash_time + dslash_space);
  node0_printf("CHECK faction: %.4g - %.4g = %.4g\n",
               faction, check, faction - check);
  dslash_space /= 3.0;
#endif

  // Check that pbg5p is purely real up to precision of inversion
  if (fabs(pbg5p.imag) > sqrt(rsqmin)) {
    node0_printf("WARNING: Im(pbg5p) = %.4g > %.4g\n",
                 pbg5p.imag, sqrt(rsqmin));
  }

  // Print results
  // Sanity check: Gaussian fermion action should be about 4N per site
//  node0_printf("FMES %.8g %.8g %.8g %.8g %.8g %.8g\n", pbp.real, pbp.imag,
//               dslash_time, dslash_space, faction, pbg5p.real);
  node0_printf("FMES %.8g %.8g %.8g %.8g\n",
               pbp.real, pbp.imag, faction, pbg5p.real);
  fflush(stdout);

  return iters;
}
// -----------------------------------------------------------------
