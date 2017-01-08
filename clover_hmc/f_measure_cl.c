// -----------------------------------------------------------------
/* Measure fermionic observables:
    psi-bar-psi, fermion action, energy and pressure, psi-bar-gamma_5-psi

   M = A - kappa*(Dslash_eo + DSLASH_oe)
   In this version, M is NOT the LU preconditioned matrix.
   MM = preconditioned matrix = A_e - kappa^2*Dslash_eo*(A_o)^{-1}*Dslash_oe

   Output is an FMES line containing the expectation value of
   the real part of Tr(1/M),
   the imaginary part of  Tr(1/M),
   the  trace of D_slash(t),
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

   psi-bar-psi = (nflavors/2) * 4 * kappa * $2
   entropy = (nflavors/2) * 2 * kappa * (-$4 + $5)
   energy = (nflavors/2) * 2 * kappa * (-$4 + (4*3 - $2)*tderiv )
   pressure = (nflavors/2) * 2 * kappa * ($5 - (4*3 - $2)*3*sderiv )
    where tderiv = partial 1/kappa_c / partial alpha_t
    where sderiv = partial 1/kappa_c / partial alpha_s

   These are the entropy, etc. summed over color and flavor.
*/
#include "cl_dyn_includes.h"

// Requires LU preconditioning in congrad
#ifndef LU
  #error "f_measure_cl requires LU preconditioning"
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int f_measure_cl() {
  register int i, j, k, dir, iters;
  register site *s;
  Real faction, dslash_time, dslash_space;
  register complex cc;
  complex pbp, pbg5p;
  half_wilson_vector hwv0,hwv1;
  wilson_vector wv0;
  msg_tag *tag0,*tag1;

  // Gaussian random vector
  FORALLSITES(i, s) {
    for(k=0;k<4;k++)for(j=0;j<DIMF;j++){
#ifdef SITERAND
      s->g_rand.d[k].c[j].real = gaussian_rand_no(&(s->site_prn));
      s->g_rand.d[k].c[j].imag = gaussian_rand_no(&(s->site_prn));
#else
      s->g_rand.d[k].c[j].real = gaussian_rand_no(&node_prn);
      s->g_rand.d[k].c[j].imag = gaussian_rand_no(&node_prn);
#endif
    }
  }

  // Copy gaussian source to chi[0]
  FORALLSITES(i,s)
    copy_wvec(&(s->g_rand), &(s->chi[0]));

  // Invert M, put result in psi[0]

  /* Load inversion control structure */
  qic.start_flag = 0;   /* Use zero initial guess for psi */

  // Load Dirac matrix parameters,
  // including temporaries specific to inverter
#ifdef BI
  qic.wv1 = F_OFFSET(tmp);
  qic.wv2 = F_OFFSET(mp);
  qic.wv3 = F_OFFSET(p);
  qic.wv4 = F_OFFSET(sss);
  iters = wilson_invert(F_OFFSET(chi[0]), F_OFFSET(psi[0]), F_OFFSET(r),
                        bicgilu_cl, &qic, (void *)&dcp);
#else
  qic.wv1 = F_OFFSET(tmp);
  qic.wv2 = F_OFFSET(mp);
  iters = wilson_invert(F_OFFSET(chi[0]), F_OFFSET(psi[0]), F_OFFSET(r),
                        cgilu_cl, &qic, (void *)&dcp);
#endif

  // DEBUG CHECK:
  /* Multiply by M and see if I get g_rand back */
  /* use dir as flag*/
  /*
     dslash_w_site(F_OFFSET(psi[0]), F_OFFSET(mp), PLUS, EVENANDODD);
     FORALLSITES(i, s)scalar_mult_add_wvec(&(s->psi[0]), &(s->mp), -kappa, &(s->mp));
     FORALLSITES(i, s){
     for(dir=0,j=0;j<4;j++)for(k=0;k<3;k++){
     if(s->g_rand.d[j].c[k].real - s->mp.d[j].c[k].real > 2e-5)dir=1;
     if(s->g_rand.d[j].c[k].imag - s->mp.d[j].c[k].imag > 2e-5)dir=1;
     if(dir)printf("%d %d %d  (%.4e , %.4e)  ( %.4e , %.4e )\n",
     i,j,k,s->g_rand.d[j].c[k].real,s->g_rand.d[j].c[k].imag,
     s->mp.d[j].c[k].real,s->mp.d[j].c[k].imag);
     }
     }
     */

  pbp = cmplx(0.0, 0.0);
  pbg5p = cmplx(0.0, 0.0);
  faction = dslash_time = dslash_space = 0.0;

  /* psi-bar-psi = g_rand.psi */
  /* psi-bar-gamma-5 psi = g_rand. gamma-5 psi */
  FORALLSITES(i, s) {
    cc = wvec_dot(&(s->g_rand), &(s->psi[0]));
    CSUM(pbp,cc);
    mult_by_gamma(&(s->psi[0]), &wv0, GAMMAFIVE);
    cc = wvec_dot(&(s->g_rand), &wv0);
    CSUM(pbg5p, cc);
  }

  /* fermion energy and pressure */
  FORALLUPDIR(dir) {
    /* multiply g_rand by one component of Dslash_adjoint, result in p.
       dot product with psi */

    /* multiply g_rand by one component of Dslash_adjoint, result in p */
    FORALLSITES(i,s){
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
    FORALLSITESDOMAIN(i,s){
      mult_su3_mat_hwvec(&(s->link[dir]),
          (half_wilson_vector *)(gen_pt[0][i]), &hwv0 );
      wp_grow(&hwv0, &(s->p), dir, MINUS);
      wp_grow((half_wilson_vector *)(gen_pt[1][i]), &wv0, dir, PLUS );
      add_wilson_vector(&wv0, &(s->p), &(s->p));
    }
    cleanup_gather(tag0);
    cleanup_gather(tag1);

    /* dot product with psi, result into energy or pressure */
    FORALLSITESDOMAIN(i,s) {
      cc = wvec_dot(&(s->psi[0]), &(s->p));
      if (dir==TUP)
        dslash_time += cc.real;
      else
        dslash_space += cc.real;
    }
  }
  g_floatsum(&dslash_time);
  g_floatsum(&dslash_space);
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
  node0_printf("FMES %.8g %.8g %.8g %.8g %.8g %.8g\n",
               (double)pbp.real, (double)pbp.imag, (double)dslash_time,
               (double)dslash_space, (double)faction, (double)pbg5p.real);
  fflush(stdout);
  return iters;
}
// -----------------------------------------------------------------
