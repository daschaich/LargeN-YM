// -----------------------------------------------------------------
#include "cl_dyn_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Update the momenta
void update_anti_hermitian(site *s, int dir, Real eps, su3_matrix_f *force) {
  su3_matrix_f tempf1, tempf2;

  uncompress_anti_hermitian(&(s->mom[dir]), &tempf1);
  scalar_mult_sub_su3_matrix_f(&tempf1, force, eps, &tempf2);
  make_anti_hermitian(&tempf2, &(s->mom[dir]));
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Update the momenta with the gauge force
double gauge_force(Real eps) {
  register int i, dir, dir2;
  register site *s;
  register Real eb3 = eps * beta / (Real)NCOL;
  double norm = 0.0;
  msg_tag *tag0, *tag1, *tag2;
  su3_matrix_f tmat1, tmat2;

  // Loop over directions, update mom[dir]
  FORALLUPDIR(dir) {
    // Initialize staple, recast to su3_matrix_f
    FORALLSITES(i, s)
      clear_su3mat_f((su3_matrix_f *)&(s->staple));

    /* Loop over other directions, computing force from plaquettes in
       the dir,dir2 plane */
    FORALLUPDIR(dir2) {
      if (dir2 == dir)
        continue;

      /* get link[dir2] from direction dir */
      tag0 = start_gather_site(F_OFFSET(linkf[dir2]),
          sizeof(su3_matrix_f), dir, EVENANDODD, gen_pt[0]);

      /* Start gather for the "upper staple" */
      tag2 = start_gather_site(F_OFFSET(linkf[dir]),
          sizeof(su3_matrix_f), dir2, EVENANDODD, gen_pt[2]);

      /* begin the computation "at the dir2DOWN point", we will
         later gather the intermediate result "to the home point" */

      /* Note: For SF we don't care if we get this wrong for
         dir<TUP and t=0, since then the staple will not be used,
         as those links and momenta are frozen */

      wait_gather(tag0);
      FORALLSITES(i, s) {
        mult_su3_an_f(&(s->linkf[dir2]), &(s->linkf[dir]), &tmat1);
        mult_su3_nn_f(&tmat1, (su3_matrix_f *)gen_pt[0][i],
                      (su3_matrix_f *)&(s->tempmat1));
      }

      /* Gather lower staple "up to home site" */
      tag1 = start_gather_site(F_OFFSET(tempmat1), sizeof(su3_matrix_f),
          OPP_DIR(dir2), EVENANDODD, gen_pt[1]);

      /*  The "upper" staple.  Note that
          one of the links has already been gathered, since it
          was used in computing the "lower" staple of the site
          above us (in dir2) */
      wait_gather(tag2);
      FORALLSITES(i, s) {
        mult_su3_nn_f(&(s->linkf[dir2]),
            (su3_matrix_f *)gen_pt[2][i], &tmat1);
        mult_su3_na_f(&tmat1, (su3_matrix_f *)gen_pt[0][i],  &tmat2);
        add_su3_matrix_f((su3_matrix_f *)&(s->staple), &tmat2,
            (su3_matrix_f *)&(s->staple));
      }

      wait_gather(tag1);

      FORALLSITES(i, s) {
        add_su3_matrix_f((su3_matrix_f *)&(s->staple),
                         (su3_matrix_f *)gen_pt[1][i],
                         (su3_matrix_f *)&(s->staple));
      }
      cleanup_gather(tag0);
      cleanup_gather(tag1);
      cleanup_gather(tag2);
    }

    // Now multiply the staple sum by the link, then update momentum 
    FORALLSITES(i, s) {
      mult_su3_na_f(&(s->linkf[dir]), (su3_matrix_f *)&(s->staple), &tmat1);
      update_anti_hermitian(s, dir, eb3, &tmat1);
      norm += (double)realtrace_su3_f(&tmat1, &tmat1);
    }
  }
  g_doublesum(&norm);
  return (eb3 * sqrt(norm) / volume);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Fermion force from hopping and clover terms
void delta(field_offset psi, field_offset p) {
  register int i, mu, nu;
  register site *s;
  msg_tag *tag0,*tag1;
  wilson_vector tvec1,tvec2;
  half_wilson_vector hvec;
  su3_matrix temp1, temp2;
  Real local_CKU0 = -clov_c * kappa / (8.0 * u0 * u0 * u0);

  FORALLUPDIR(mu) {
#ifdef CG_DEBUG
    node0_printf("mu = %d\n",mu);
#endif

    FORALLSITES(i, s) {
      clear_su3mat(&(s->tempmat1));
      wp_shrink((wilson_vector *)F_PT(s, psi), &(s->htmp[0]), mu, PLUS);
      wp_shrink((wilson_vector *)F_PT(s, p), &(s->htmp[1]), mu, MINUS);
    }
    tag0 = start_gather_site(F_OFFSET(htmp[0]), sizeof(half_wilson_vector),
        mu, EVENANDODD, gen_pt[0]);
    tag1 = start_gather_site(F_OFFSET(htmp[1]), sizeof(half_wilson_vector),
        mu, EVENANDODD, gen_pt[1]);
    wait_gather(tag0);
    wait_gather(tag1);

    /* First the U d(dslash)/dU terms (do for only one value of nu) */
    FORALLSITES(i, s) {
      /* psi and p parallel transported in from positive directions */
      mult_su3_mat_hwvec(&(s->link[mu]),
                         (half_wilson_vector *)gen_pt[0][i], &hvec);
      wp_grow(&hvec, &tvec1, mu, PLUS);
      /* i even => tvec1 = (1+gamma_mu)*U*Aodd^(-1)*D*psi,
         i odd  => tvec1 = (1+gamma_mu)*U*psi */

      mult_su3_mat_hwvec(&(s->link[mu]),
          (half_wilson_vector *)gen_pt[1][i], &hvec);
      wp_grow(&hvec, &tvec2, mu, MINUS);
      /* i even => tvec2 = (1-gamma_mu)*U*Aodd^(-1)*D_adj*M*psi,
         i odd  => tvec2 = (1-gamma_mu)*U*M*psi */

      su3_projector_w(&tvec1, (wilson_vector *)F_PT(s,p), &temp1);
      su3_projector_w(&tvec2, (wilson_vector *)F_PT(s,psi), &temp2);
      add_su3_matrix(&temp1, &temp2, &(s->tempmat1));
      scalar_mult_su3_matrix(&(s->tempmat1), -kappa, &(s->tempmat1));
    }
    cleanup_gather(tag0);
    cleanup_gather(tag1);

    // Now the U dA/dU terms
    FORALLUPDIR(nu) {
      if (nu == mu)
        continue;

      // U dA/dU from U dM/dU
      udadu_mu_nu(p, psi, F_OFFSET(tempmat2), mu, nu, EVENANDODD);
      FORALLSITES(i, s) {
        scalar_mult_add_su3_matrix(&(s->tempmat1), &(s->tempmat2),
                                   local_CKU0, &(s->tempmat1));
      }

      // U dA/dU from U dM^dagger/dU
      udadu_mu_nu(psi,p, F_OFFSET(tempmat2), mu, nu, EVENANDODD);
      FORALLSITES(i, s) {
        scalar_mult_add_su3_matrix(&(s->tempmat1), &(s->tempmat2),
                                   local_CKU0, &(s->tempmat1));
      }
    }
    FORALLSITES(i, s)
      add_su3_matrix(&(s->Force[mu]), &(s->tempmat1), &(s->Force[mu]));
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Gauge force from fermion irrep plaquette term
void gauge_force_frep(field_offset scratch1, field_offset scratch2,
                      field_offset forcemat, int dir1) {

  register int i, dir2;
  register site *s;
  msg_tag *tag0,*tag1, *tag2;
  su3_matrix tmat1,tmat2;
  /* divide by nflavors=2 to compensate for
     ferm_epsilon = nflavors*eps = 2*eps
     */
  Real bfrep_force = beta_frep/(Real)DIMF/(Real)nflavors;

  /* Loop over directions for updating mom[dir1] is in the calling routine.
     Here, we loop over other directions, computing force from frep plaquettes
     in the dir1,dir2 plane
     */
  /* initialize staple */
  FORALLSITES(i, s)
    clear_su3mat((su3_matrix *)F_PT(s, scratch1));

  for (dir2=XUP;dir2<=TUP;dir2++)if (dir2 != dir1) {
    /* get link[dir2] from direction dir1 */
    tag0 = start_gather_site(F_OFFSET(link[dir2]),
        sizeof(su3_matrix), dir1, EVENANDODD, gen_pt[0]);

    /* Start gather for the "upper staple" */
    tag2 = start_gather_site(F_OFFSET(link[dir1]),
        sizeof(su3_matrix), dir2, EVENANDODD, gen_pt[2]);

    /* begin the computation "at the dir2DOWN point", we will
       later gather the intermediate result "to the home point" */

    wait_gather(tag0);
    FORALLSITES(i, s) {
      mult_su3_an(&(s->link[dir2]), &(s->link[dir1]), &tmat1);
      mult_su3_nn(&tmat1, (su3_matrix *)gen_pt[0][i],
          ((su3_matrix *)F_PT(s,scratch2)));
    }

    /* Gather lower staple "up to home site" */
    tag1 = start_gather_site(scratch2, sizeof(su3_matrix),
        OPP_DIR(dir2), EVENANDODD, gen_pt[1]);

    /*  The "upper" staple.  Note that
        one of the links has already been gathered, since it
        was used in computing the "lower" staple of the site
        above us (in dir2) */
    wait_gather(tag2);
    FORALLSITES(i, s) {
      mult_su3_nn(&(s->link[dir2]),
          (su3_matrix *)gen_pt[2][i], &tmat1);
      mult_su3_na(&tmat1, (su3_matrix *)gen_pt[0][i],  &tmat2);
      add_su3_matrix(((su3_matrix *)F_PT(s,scratch1)), &tmat2,
          ((su3_matrix *)F_PT(s,scratch1)));
    }

    wait_gather(tag1);

    FORALLSITES(i, s) {
      add_su3_matrix((su3_matrix *)F_PT(s, scratch1),
                     (su3_matrix *)gen_pt[1][i],
                     (su3_matrix *)F_PT(s, scratch1));
    }
    cleanup_gather(tag0);
    cleanup_gather(tag1);
    cleanup_gather(tag2);
  }

  // Now multiply the staple sum by the link, then add to force
  FORALLSITES(i, s) {
    mult_su3_na(&(s->link[dir1]), (su3_matrix *)F_PT(s, scratch1), &tmat1);
    scalar_mult_add_su3_matrix((su3_matrix *)F_PT(s, forcemat), &tmat1,
                               bfrep_force, (su3_matrix *)F_PT(s, forcemat));
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Update the  momenta with the fermion force
/* Assumes that the conjugate gradient has been run, with the answer in psi */
/* pseudofermion vector is chi */

/* Use the LU preconditioned fermion matrix, where
  the fermion spinors live on even sites only.  In other words, if
  Dslash_oe is the dslash operator with its source on even sites and
  its result on odd sites, etc.:
    M = A_e - kappa^2 * Dslash_eo (A_o)^{-1}* Dslash_oe
*/
double fermion_force(Real eps1, Real eps2) {
  register int i, mu, nu;
  register site *s;
  int level;
  Real local_CKU0 = -clov_c * kappa / (8.0 * u0 * u0 * u0);
  Real ferm_epsilon = nflavors * eps1, tr;
  double tmpnorm, maxnorm = 0.0, norm = 0.0, returnit = 0.0;
  su3_matrix temp1;
  su3_matrix_f tempf1;

#ifdef TIMING
TIC(1)
#endif

  // Clear the force collectors
  FORALLUPDIR(mu) {
    FORALLSITES(i, s)
      clear_su3mat(&(s->Force[mu]));
  }

  /* while we're at it, get the diagonal clover entry */
  /* With A_odd == sigma_mu_nu F_mu_nu, now add force term
     U dF_mu_nu/dU Tr_{dirac} (sigma_mu_nu A_odd^{-1}) */
  FORALLUPDIR(mu) {
    FORALLUPDIR(nu) {
      if (nu == mu)
        continue;
      // Use staple as a temporary su3_matrix
      tr_sigma_ldu_mu_nu_site(F_OFFSET(staple), mu, nu);
      udadu_mat_mu_nu(F_OFFSET(staple), F_OFFSET(tempmat2), mu, nu);
      FORALLSITES(i, s) {
        scalar_mult_add_su3_matrix(&(s->Force[mu]), &(s->tempmat2),
                                   2.0 * local_CKU0, &(s->Force[mu]));
      }
    }
  }

  FORALLUPDIR(mu) {
    gauge_force_frep(F_OFFSET(staple), F_OFFSET(tempmat2),
                     F_OFFSET(Force[mu]), mu);
  }

    if (num_masses == 1)
      level = 0;
    else
      level = 1;
    prepare_vecs(level);
    delta(F_OFFSET(psi[level]), F_OFFSET(p));

    if (eps2 > 0.0 && num_masses == 2) {
      dslash_w_site(F_OFFSET(psi[0]), F_OFFSET(tmp), PLUS, ODD);
      mult_ldu_site(F_OFFSET(tmp), F_OFFSET(psi[0]), ODD);
      FORODDSITES(i, s)
        scalar_mult_wvec(&(s->psi[0]), kappa, &(s->psi[0]));

      dslash_w_site(F_OFFSET(psi[0]), F_OFFSET(p), PLUS, EVEN);
      mult_ldu_site(F_OFFSET(psi[0]), F_OFFSET(tmp), EVEN);
      FOREVENSITES(i, s)
        scalar_mult_add_wvec(&(s->tmp), &(s->p), -kappa, &(s->p));

      dslash_w_site(F_OFFSET(p), F_OFFSET(tmp), MINUS, ODD);
      mult_ldu_site(F_OFFSET(tmp), F_OFFSET(p), ODD);
      FORODDSITES(i, s)
        scalar_mult_wvec(&(s->p), kappa, &(s->p));

      /* M = A_even - kappa^2 * Dslash * A_odd^{-1} * Dslash
   psi(even) = psi(even)
   psi(odd)  = kappa * A_odd^{-1} * Dslash * psi(even)
   p(even) = M * psi(even)
   p(odd)  = kappa * A_odd^{-1} * Dslash_adjoint * M * psi(even). */

      /* And the overall factor of shift^2, here we also adjust the stepsize */
      tr = shift * shift * eps2 / eps1;
      FORALLSITES(i,s)
        scalar_mult_wvec(&(s->p), tr, &(s->p));

      delta(F_OFFSET(psi[0]), F_OFFSET(p));
    }

      FORALLUPDIR(mu) {
        FORALLSITES(i, s) {
          /* back from fermion irrep to fundamental links
             compatible with fat links
             W = link in fermion irrep (can be thin or fat)
             V = fat fundamental link
             U = thin fundamental link
             */
          /* remove W from W*dH/dW^T      */
          mult_su3_an(&(s->link[mu]), &(s->Force[mu]), &temp1);

          /* dH/dV^T from dW/dV^T dH/dW^T */
          chain_rule(&tempf1, &temp1, gauge_field[mu] + i);

          // Take care of boundary conditions
          apply_bc(&tempf1, mu, s->t);

          // Not done yet.  For now, save dH/dV^T
          su3mat_copy_f(&tempf1, Sigma[mu]+i);
        }
      }

    // More chain rule: dH/dU^T = dV/dU^T dH/dV^T
    // stout_force receives and returns force in Sigma[mu]
    stout_force1();
    // Finally we can update
    FORALLUPDIR(mu) {
      FORALLSITES(i, s) {
        // First multiply back U*dH/dU^T
        mult_su3_nn_f(&(s->linkf[mu]), Sigma[mu]+i, &tempf1);
        // Now update
        update_anti_hermitian(s, mu, ferm_epsilon, &tempf1);
        tmpnorm = (double)realtrace_su3_f(&tempf1, &tempf1);
        norm += tmpnorm;
        if (tmpnorm > maxnorm)
          maxnorm = tmpnorm;
      }
    }
      g_doublesum(&norm);
      returnit=ferm_epsilon*sqrt(norm)/volume*2;

#ifdef TIMING
  TOC(1, time_fermion_force)
#endif

/**dtime += dclock();
if (this_node==0)printf("F_FORCE: time = %e mflops = %e\n",
dtime, (double)(5584.0*volume/(1.0e6*dtime*numnodes())));**/
    return(returnit);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void prepare_vecs(int level) {
  register int i;
  register site *s;
  complex tc;
  wilson_vector tvec, tvec2;

  if (level != 0)
    tc = cmplx(0, shift);

  dslash_w_site(F_OFFSET(psi[level]), F_OFFSET(tmp), PLUS, ODD);
  mult_ldu_site(F_OFFSET(tmp), F_OFFSET(psi[level]), ODD);
  FORODDSITES(i, s)
    scalar_mult_wvec(&(s->psi[level]), kappa, &(s->psi[level]));

  dslash_w_site(F_OFFSET(psi[level]), F_OFFSET(p), PLUS, EVEN);
  mult_ldu_site(F_OFFSET(psi[level]), F_OFFSET(tmp), EVEN);
  FOREVENSITES(i, s) {
    if (level == 1) {
      scalar_mult_add_wvec(&(s->tmp), &(s->p), -kappa, &tvec);
      /* that was M*psi now we need to add i*shift*gamma_5*p */
      mult_by_gamma(&(s->psi[level]), &tvec2, GAMMAFIVE);
      c_scalar_mult_add_wvec(&tvec, &tvec2, &tc, &(s->p));
    }
    else
      scalar_mult_add_wvec(&(s->tmp), &(s->p), -kappa, &(s->p));
  }
  dslash_w_site(F_OFFSET(p), F_OFFSET(tmp), MINUS, ODD);
  mult_ldu_site(F_OFFSET(tmp), F_OFFSET(p), ODD);
  FORODDSITES(i, s)
    scalar_mult_wvec(&(s->p), kappa, &(s->p));

  /* M = A_even - kappa^2 * Dslash * A_odd^{-1} * Dslash
     psi(even) = psi(even)
     psi(odd)  = kappa * A_odd^{-1} * Dslash * psi(even)
     p(even) = M * psi(even)
     p(odd)  = kappa * A_odd^{-1} * Dslash_adjoint * M * psi(even). */
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void update_h(Real eps) {
  gauge_force(eps);
  fermion_force(eps, 0.0);
}
// -----------------------------------------------------------------
