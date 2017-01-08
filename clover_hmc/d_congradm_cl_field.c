// Solve Mdag.M psi = chi for clover fermions for dynamical HMC updates

/* Use the LU preconditioned fermion matrix, where
  the fermion spinors live on even sites only.  In other words, if
  Dslash_oe is the dslash operator with its source on even sites and
  its result on odd sites, etc.:
   M = A_e - kappa^2 * Dslash_eo * (A_o)^{-1} * Dslash_oe
*/

/* adapted to "field" instead of "site" */

/* This version looks at the initial vector every "niter" passes */
/* The source vector is in "chi", and the initial guess and answer
   in "psi".  "r" is the residual vector, and "p" and "mp" are
   working vectors for the conjugate gradient.
   niter = maximum number of iterations.
   rsqmin = desired rsq, quit when we reach rsq = rsqmin*source_norm.
*/

#include "cl_dyn_includes.h"

int congrad_cl_m(int niter, Real rsqmin, Real *final_rsq_ptr,
                 field_offset src, field_offset dest, Real mshift) {

  register int i;
  register site *s;
  int iteration;
  Real a, b;
  double source_norm, rsqstop;
  double rsq,oldrsq,pkp;	/* Sugar's a,b,resid**2,previous resid*2 */
  /* pkp = cg_p.K.cg_p */
  void dslash_w_field_special();
  msg_tag *tag[8], *tag2[8];
  Real kappaSq = -kappa * kappa, MSq = mshift * mshift;

#ifdef TIMING
  TIC(0)
#endif

    /* create space for fields, copy source and initial guess */
    wilson_vector *psi, *chi, *tmp, *mp, *p, *r;
    FIELD_ALLOC(psi,wilson_vector)
    FIELD_ALLOC(chi,wilson_vector)
    FIELD_ALLOC(tmp,wilson_vector)
    FIELD_ALLOC(mp,wilson_vector)
    FIELD_ALLOC(p,wilson_vector)
    FIELD_ALLOC(r,wilson_vector)

    FORALLSITES(i,s) {
      copy_wvec((wilson_vector *)F_PT(s,src),&(chi[i]));
      copy_wvec((wilson_vector *)F_PT(s,dest),&(psi[i]));
    }

  double dtime;
  dtime= -dclock();

  iteration=0;
start:
  /* mp <-  M_adjoint*M*psi
     r,p <- chi - mp
     rsq = |r|^2
     source_norm = |chi|^2
     */
  rsq = source_norm = 0.0;
  mult_ldu_field(psi, tmp, EVEN);
  dslash_w_field_special(psi, mp, PLUS, ODD, tag, 0);
  mult_ldu_field(mp, tmp, ODD);
  dslash_w_field_special(tmp, mp, PLUS, EVEN, tag2, 0);
  FOREVENSITESDOMAIN(i,s) {
    scalar_mult_add_wvec(tmp+i, mp+i, kappaSq, mp+i);
  }/* So far this was the application of M */

  mult_ldu_field(mp, tmp, EVEN);
  dslash_w_field_special(mp, mp, MINUS, ODD, tag, 1);
  mult_ldu_field(mp, tmp, ODD);
  dslash_w_field_special(tmp, mp, MINUS, EVEN, tag2, 1);
  FOREVENSITESDOMAIN(i,s) {
    scalar_mult_add_wvec(tmp+i, mp+i, kappaSq, mp+i);
    /* Now we have Mdagger*M applied to psi and have the result in mp
       in the next step we add the shift part, i.e. mshift^2*psi */
    /*  add mshift^2 */
    if (mshift > 0.0)
      scalar_mult_add_wvec(mp + i, psi + i, MSq, mp + i);

    sub_wilson_vector(chi + i, mp + i, r + i);
    p[i] = r[i];
    source_norm += (double)magsq_wvec(chi + i);
    rsq += (double)magsq_wvec(r + i);
  }
  g_doublesum(&source_norm);
  g_doublesum(&rsq);
  iteration++ ;	/* iteration counts number of multiplications
                   by M_adjoint*M */
  total_iters++;
#ifdef CG_DEBUG
  node0_printf("CG source_norm = %.4g\n", source_norm);
  node0_printf("CG iter %d, rsq %.4g, pkp %.4g, a %.4g\n",
               iteration, (double)rsq, (double)pkp, (double)a);
#endif

  rsqstop = rsqmin * source_norm;
  if (rsq <= rsqstop) {
    *final_rsq_ptr= (Real)rsq;
    FORALLUPDIR(i) {
      cleanup_gather(tag[i]);
      cleanup_gather(tag2[i]);
      cleanup_gather(tag[OPP_DIR(i)]);
      cleanup_gather(tag2[OPP_DIR(i)]);
    }

    cleanup_dslash_temps();
    cleanup_tmp_links();

    FORALLSITES(i, s)
      copy_wvec(&(psi[i]), (wilson_vector *)F_PT(s,dest));

    free(psi);
    free(chi);
    free(tmp);
    free(mp);
    free(p);
    free(r);

#ifdef TIMING
    TOC(0,time_dcongrad)
#endif

      return (iteration);
  }

  /* main loop - do until convergence or time to restart */
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
    FOREVENSITESDOMAIN(i,s)
      scalar_mult_add_wvec(tmp + i, mp + i, kappaSq, mp + i);

    mult_ldu_field(mp, tmp, EVEN);
    dslash_w_field_special(mp, mp, MINUS, ODD, tag, 1);
    mult_ldu_field(mp, tmp, ODD);
    dslash_w_field_special(tmp, mp, MINUS, EVEN, tag2, 1);
    FOREVENSITESDOMAIN(i,s) {
      scalar_mult_add_wvec(tmp+i, mp + i, kappaSq, mp + i);
      /*  add mshift^2 */
      if (mshift > 0.0)
        scalar_mult_add_wvec(mp + i, p + i, MSq, mp + i);
      pkp += (double)wvec_rdot(p + i, mp + i);
    }
    g_doublesum(&pkp);
    iteration++;
    total_iters++;

    a = (Real)(rsq / pkp);
    rsq = 0.0;
    FOREVENSITESDOMAIN(i, s) {
      scalar_mult_add_wvec(psi + i, p + i, a, psi + i);
      scalar_mult_add_wvec(r + i, mp + i, -a, r + i);
      rsq += (double)magsq_wvec(r + i);
    }
    g_doublesum(&rsq);
    /* if (this_node==0)printf("congrad2: iter %d, rsq %e, pkp %e, a %e\n",
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
      /* if (this_node==0)printf("CONGRAD2: time = %e iters = %d mflops = %e\n",
         dtime,iteration,(double)(2840.0*volume*iteration/(1.0e6*dtime*numnodes()))); */
      cleanup_dslash_temps();
      cleanup_tmp_links();

      FORALLSITES(i,s) {
        copy_wvec(&(psi[i]),(wilson_vector *)F_PT(s,dest));
      }
      free(psi); free(chi); free(tmp);
      free(mp); free(p); free(r);

#ifdef TIMING
      TOC(0,time_dcongrad)
#endif

      return iteration;
    }

    b = (Real)(rsq/oldrsq);
    FOREVENSITESDOMAIN(i,s)
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
