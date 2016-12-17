/* SF:
The SF version is based on
clover_dynamical/update_h_cl.c  for the fermion force, and
generic_pg/update_h.c           for the gauge force.
*/

/*
Note: If the plaquette action is used,
this routine is valid also for clover fermions (with the standard (a.)p.b.c.)
and can replace update_h_cl.c .  See also d_action.c .
*/
#include "cl_dyn_includes.h"
double gauge_force(Real eps);
#ifdef BETA_FREP
void gauge_force_frep(field_offset scratch1, field_offset scratch2,
                      field_offset forcemat, int dir1);
#endif
void update_anti_hermitian( site *st, int dir, Real eps, su3_matrix_f *force);

/*** update_h ***************************************************************/
void update_h(Real eps){
  double normg,normf;
  normg=normf=0.0;
    /* gauge field force */
    normg +=gauge_force(eps);
    /* fermionic force */
#ifdef NHYP
    sf_coupling_flag=FORCE;
#endif
    normf += fermion_force(eps,0.0);
} /* update_h */


/* update the momenta with the gauge force */

/* SF gauge force.  Based on generic_pg/update_h.c.
Originally adapted for SF by UMH, Jan 2000
Here adapted to new format of SF boundary values, and to NCOL. */

/*** gauge-action force *****************************************************/
double gauge_force(Real eps) {
register int i,dir1,dir2;
register site *st;
msg_tag *tag0,*tag1,*tag2;
su3_matrix_f tmat1,tmat2;
register Real eb3;
#ifdef VERBOSE
double tmpnorm,maxnorm=0.0;
#endif
 double norm=0.0;


    eb3 = eps*beta/(Real)NCOL;
    /* Loop over directions, update mom[dir1] */



    for(dir1=XUP; dir1<=TUP; dir1++){
	/* Loop over other directions, computing force from plaquettes in
	   the dir1,dir2 plane */
	FORALLSITES(i,st){
        /* initialize staple, recast to su3_matrix_f */
            clear_su3mat_f( (su3_matrix_f *)&(st->staple) );
        }
	for(dir2=XUP;dir2<=TUP;dir2++)if(dir2 != dir1){

	    /* get link[dir2] from direction dir1 */
	    tag0 = start_gather_site( F_OFFSET(linkf[dir2]),
                sizeof(su3_matrix_f), dir1, EVENANDODD, gen_pt[0] );

	    /* Start gather for the "upper staple" */
	    tag2 = start_gather_site( F_OFFSET(linkf[dir1]),
                sizeof(su3_matrix_f), dir2, EVENANDODD, gen_pt[2] );

	    /* begin the computation "at the dir2DOWN point", we will
		later gather the intermediate result "to the home point" */

	    /* Note: For SF we don't care if we get this wrong for
		     dir1<TUP and t=0, since then the staple will not be used,
		     as those links and momenta are frozen */

	    wait_gather(tag0);
            /* SF: supersede the gather with boundary links if necessary */
            FORALLSITES(i,st) {
                gen_pt[0][i] = CHOOSE_NBR(i,st,dir1,linkf_bndr_up[dir2],0);
            }

	    FORALLSITES(i,st){
		mult_su3_an_f( &(st->linkf[dir2]), &(st->linkf[dir1]), &tmat1);
		mult_su3_nn_f( &tmat1, (su3_matrix_f *)gen_pt[0][i],
		   (su3_matrix_f *)&(st->tempmat1) );
	    }

	    /* Gather lower staple "up to home site" */
	    tag1 = start_gather_site( F_OFFSET(tempmat1), sizeof(su3_matrix_f),
		OPP_DIR(dir2), EVENANDODD, gen_pt[1] );

	    /*  The "upper" staple.  Note that
		one of the links has already been gathered, since it
		was used in computing the "lower" staple of the site
		above us (in dir2) */
	    wait_gather(tag2);
            /* SF: supersede the gather with boundary links if necessary */
            FORALLSITES(i,st) {
                gen_pt[2][i] = CHOOSE_NBR(i,st,dir2,linkf_bndr_up[dir1],2);
            }

	    FORALLSITES(i,st){
	        mult_su3_nn_f( &(st->linkf[dir2]),
       		   (su3_matrix_f *)gen_pt[2][i], &tmat1);
		mult_su3_na_f( &tmat1, (su3_matrix_f *)gen_pt[0][i],  &tmat2);
		add_su3_matrix_f( (su3_matrix_f *)&(st->staple), &tmat2,
                   (su3_matrix_f *)&(st->staple) );
	    }

	    wait_gather(tag1);

	    FORALLDYNLINKS(i,st,dir1){
		add_su3_matrix_f( (su3_matrix_f *)&(st->staple),
                  (su3_matrix_f *)gen_pt[1][i], (su3_matrix_f *)&(st->staple));
	    }
	    cleanup_gather(tag0);
	    cleanup_gather(tag1);
	    cleanup_gather(tag2);
	}
	/* Now multiply the staple sum by the link, then update momentum */

	FORALLDYNLINKS(i,st,dir1){
	    mult_su3_na_f( &(st->linkf[dir1]),
                (su3_matrix_f *)&(st->staple), &tmat1 );
            update_anti_hermitian(st, dir1, eb3, &tmat1);

	    norm += (double)realtrace_su3_f( &tmat1 , &tmat1);
#ifdef VERBOSE
        /* If we desired, also look for the maximal and minimal force */
        tmpnorm = (double)realtrace_su3_f( &tmat1 , &tmat1);
        if (tmpnorm>maxnorm) {maxnorm=tmpnorm;}
#endif



	}
    }

#ifdef VERBOSE
g_doublemax(&maxnorm);
node0_printf("GAUGEFORCEMAX %e\n",eb3*sqrt(maxnorm));
#endif

g_doublesum(&norm);
return (eb3*sqrt(norm)/volume);

} /* gauge_force */



/* update the  momenta with the fermion force */
/* Assumes that the conjugate gradient has been run, with the answer in psi */
/* pseudofermion vector is chi */

/* if "LU" is defined use the LU preconditioned fermion matrix, where
  the fermion spinors live on even sites only.  In other words, if
  Dslash_oe is the dslash operator with its source on even sites and
  its result on odd sites, etc.:

  without LU:
   	M = A - kappa Dslash
  with LU:
  	M = A_e - kappa^2 * Dslash_eo (A_o)^{-1}* Dslash_oe

*/

/*** fermion-det force *****************************************************/
double fermion_force(Real eps1, Real eps2){
  register int i, mu,nu;
  register site *st;
  su3_matrix temp1;
  su3_matrix_f tempf1;
  int level;
  Real local_CKU0=-clov_c*kappa/(8.0*u0*u0*u0);

  void delta(field_offset psi, field_offset p);

#ifdef TIMING
TIC(1)
#endif

#ifndef NHYP
  su3_matrix_f tempf2;
#endif
  Real ferm_epsilon;
#ifdef NHYP
  double tmpnorm,maxnorm=0.0;
  double norm=0;
#endif
  double returnit=0.0;

  /**double dtime,dclock();
     dtime = -dclock();**/
  ferm_epsilon = nflavors*eps1;
  /* First we have to zero the force collectors */
  for (mu = XUP; mu <= TUP; mu++) FORALLSITES(i,st) clear_su3mat(&(st->Force[mu]));

  /* while we're at it, get the diagonal clover entry */
  /* With A_odd == sigma_mu_nu F_mu_nu , now add force term
     U dF_mu_nu/dU Tr_{dirac} (sigma_mu_nu A_odd^{-1}) */
#ifdef LU
  for (mu = XUP; mu <= TUP; mu++)for(nu=XUP;nu<=TUP;nu++) if(nu!=mu){
    /* staple is used as a temporary su3_matrix */
    tr_sigma_ldu_mu_nu_site( F_OFFSET(staple), mu, nu );
    /* SF: matsrc provided to udadu_mat_mu_nu must be zero on t=0 sites */
#ifdef SF
    FORALLSITES(i,st) {
      if (st->t==0) {
	clear_su3mat( &(st->staple) );
      }
    }
#endif
    udadu_mat_mu_nu( F_OFFSET(staple), F_OFFSET(tempmat2), mu, nu );
    FORALLSITES(i,st) {
      scalar_mult_add_su3_matrix( &(st->Force[mu]),
				  &(st->tempmat2), 2.0*local_CKU0, &(st->Force[mu]) );
    }
  }
#endif /*LU*/

#ifdef BETA_FREP
  for (mu = XUP; mu <= TUP; mu++){
      gauge_force_frep(F_OFFSET(staple), F_OFFSET(tempmat2),
                         F_OFFSET(Force[mu]), mu);
  }
#endif

    if(num_masses==1) level=0;
    else level=1;
    prepare_vecs(level);
    delta(F_OFFSET(psi[level]),F_OFFSET(p));

    if(eps2 >0.0 && num_masses==2){
      /* SF: all t=0 wilson vectors stay zero */
#ifdef LU
      dslash_w_site( F_OFFSET(psi[0]), F_OFFSET(tmp), PLUS, ODD );
      mult_ldu_site( F_OFFSET(tmp), F_OFFSET(psi[0]), ODD );
      FORODDSITESDOMAIN(i,st) {
	scalar_mult_wvec( &(st->psi[0]), kappa, &(st->psi[0]) );
      }
      dslash_w_site( F_OFFSET(psi[0]), F_OFFSET(p), PLUS, EVEN );
      mult_ldu_site( F_OFFSET(psi[0]), F_OFFSET(tmp), EVEN );
      FOREVENSITESDOMAIN(i,st) {
	scalar_mult_add_wvec( &(st->tmp), &(st->p), -kappa, &(st->p) );
      }
      dslash_w_site( F_OFFSET(p), F_OFFSET(tmp), MINUS, ODD );
      mult_ldu_site( F_OFFSET(tmp), F_OFFSET(p), ODD );
      FORODDSITESDOMAIN(i,st) {
	scalar_mult_wvec( &(st->p), kappa, &(st->p) );
      }

      /* M = A_even - kappa^2 * Dslash * A_odd^{-1} * Dslash
	 psi(even) = psi(even)
	 psi(odd)  = kappa * A_odd^{-1} * Dslash * psi(even)
	 p(even) = M * psi(even)
	 p(odd)  = kappa * A_odd^{-1} * Dslash_adjoint * M * psi(even). */
#else /* no LU */
      dslash_w_site( F_OFFSET(psi[0]), F_OFFSET(p), PLUS, EVENANDODD );
      mult_ldu_site( F_OFFSET(psi[0]), F_OFFSET(tmp), EVENANDODD );
      FORALLSITESDOMAIN(i,st) {
	scalar_mult_add_wvec( &(st->tmp), &(st->p), -kappa, &(st->p) );
      }
      /* M = A - kappa * Dslash
	 psi = psi
	 p = M * psi          */
#endif /* LU */
      /* And the overall factor of shift^2, here we also adjust the stepsize */
      FORALLSITESDOMAIN(i,st) {
	scalar_mult_wvec( &(st->p), shift*shift*eps2/eps1, &(st->p) );
      }


      delta(F_OFFSET(psi[0]),F_OFFSET(p));

    } /*end eps2 part */



    for(mu=XUP;mu<=TUP;mu++){

      FORALLDYNLINKS(i,st,mu) {

	/* back from fermion irrep to fundamental links
	   compatible with fat links
	   W = link in fermion irrep (can be thin or fat)
	   V = fat fundamental link
	   U = thin fundamental link
	*/
	/* remove W from W*dH/dW^T      */
	mult_su3_an(&(st->link[mu]), &(st->Force[mu]), &temp1 );

#ifndef NHYP
	/* dH/dU^T from dW/dU^T dH/dW^T */
	chain_rule(&tempf1, &temp1, &(st->linkf[mu]) );
#else
	/* dH/dV^T from dW/dV^T dH/dW^T */
	chain_rule(&tempf1, &temp1, gauge_field[mu]+i );
#endif

	/* take care of boundary conditions */
	apply_bc(&tempf1, mu, st->t );

	/* if we have thin links, we may instead use:
	   fund_rep_force( &tempf2, &(st->Force[mu]));
	*/
#ifndef NHYP  /* thin links                 */
	/* multiply back U*dH/dU^T      */
	mult_su3_nn_f(&(st->linkf[mu]), &tempf1, &tempf2 );
	/* now update                   */
	update_anti_hermitian(st, mu, ferm_epsilon, &tempf2 );
#else
	/* NHYP: we're not done yet.  for now, save dH/dV^T
	 */
	su3mat_copy_f(&tempf1, Sigma[mu]+i);
#endif
      }
    } /* end loop over mu */

#ifdef NHYP
    /* more chain rule: dH/dU^T = dV/dU^T dH/dV^T
       stout_force receives and returns force in Sigma[mu]
    */

    stout_force1();


    switch(sf_coupling_flag){

    case FORCE:
      /* finally we can update. */
      for(mu=XUP;mu<=TUP;mu++) {
	FORALLDYNLINKS(i,st,mu) {

	  /* first multiply back U*dH/dU^T      */
	  mult_su3_nn_f(&(st->linkf[mu]), Sigma[mu]+i, &tempf1 );
	  /* now update                         */
	  update_anti_hermitian(st, mu, ferm_epsilon, &tempf1 );
	  norm += (double)realtrace_su3_f(&tempf1,&tempf1);

	  tmpnorm = (double)realtrace_su3_f(&tempf1,&tempf1);
#ifdef VERBOSE
	  node0_printf("XXXFFF %e\n",tmpnorm);
#endif
	  if (tmpnorm>maxnorm) {maxnorm=tmpnorm;}


	}
      }

#ifdef VERBOSE
      g_doublemax(&maxnorm);
      node0_printf("FERMIONFORCEMAX %e\n",ferm_epsilon*sqrt(maxnorm));
      fflush(stdout);
      exit(1);
#endif
      g_doublesum(&norm);
      returnit=ferm_epsilon*sqrt(norm)/volume*2;

      break;

    case SF_COUPLING:
      /* we're done */
      returnit=0.0;
      break;

    default:
      if(this_node==0)printf("Illegal value of sf_coupling_flag\n");
      terminate(1);
    }
#endif /* NHYP */

#ifdef TIMING
TOC(1,time_fermion_force)
#endif

/**dtime += dclock();
if(this_node==0)printf("F_FORCE: time = %e mflops = %e\n",
dtime, (double)(5584.0*volume/(1.0e6*dtime*numnodes())) );**/
    return(returnit);
} /* end fermion_force */


void prepare_vecs(int level){
  register int i;
  register site *st;
  complex ctmp;
  wilson_vector tmpvec,tmpvec2;

  if(level!=0) {ctmp=cmplx(0,shift);}

  /* SF: all t=0 wilson vectors stay zero */
#ifdef LU
  dslash_w_site( F_OFFSET(psi[level]), F_OFFSET(tmp), PLUS, ODD );
  mult_ldu_site( F_OFFSET(tmp), F_OFFSET(psi[level]), ODD );
  FORODDSITESDOMAIN(i,st) {
    scalar_mult_wvec( &(st->psi[level]), kappa, &(st->psi[level]) );
  }
  dslash_w_site( F_OFFSET(psi[level]), F_OFFSET(p), PLUS, EVEN );
  mult_ldu_site( F_OFFSET(psi[level]), F_OFFSET(tmp), EVEN );
  FOREVENSITESDOMAIN(i,st) {
    if(level==1){
      scalar_mult_add_wvec( &(st->tmp), &(st->p), -kappa, &tmpvec );
      /* that was M*psi now we need to add i*shift*gamma_5*p */
      mult_by_gamma(&(st->psi[level]),&tmpvec2,GAMMAFIVE);
      c_scalar_mult_add_wvec2(&tmpvec,&tmpvec2,ctmp,&(st->p));
    }
    else
      scalar_mult_add_wvec( &(st->tmp), &(st->p), -kappa, &(st->p) );
  }
  dslash_w_site( F_OFFSET(p), F_OFFSET(tmp), MINUS, ODD );
  mult_ldu_site( F_OFFSET(tmp), F_OFFSET(p), ODD );
  FORODDSITESDOMAIN(i,st) {
    scalar_mult_wvec( &(st->p), kappa, &(st->p) );
  }

  /* M = A_even - kappa^2 * Dslash * A_odd^{-1} * Dslash
     psi(even) = psi(even)
     psi(odd)  = kappa * A_odd^{-1} * Dslash * psi(even)
     p(even) = M * psi(even)
     p(odd)  = kappa * A_odd^{-1} * Dslash_adjoint * M * psi(even). */
#else /* no LU */
  dslash_w_site( F_OFFSET(psi[level]), F_OFFSET(p), PLUS, EVENANDODD );
  mult_ldu_site( F_OFFSET(psi[level]), F_OFFSET(tmp), EVENANDODD );
  FORALLSITESDOMAIN(i,st) {
    if(level==1){
      scalar_mult_add_wvec( &(st->tmp), &(st->p), -kappa, &tmpvec );
      /* that was M*psi now we need to add i*shift*gamma_5*p */
      mult_by_gamma(&(s->psi[level]),&tmpvec2,GAMMAFIVE);
      c_scalar_mult_add_wvec2(&tmpvec,&tmpvec2,ctmp,&(st->p));
    }
    else
      scalar_mult_add_wvec( &(st->tmp), &(st->p), -kappa, &(st->p) );

  }
  /* M = A - kappa * Dslash
     psi = psi
     p = M * psi          */
#endif /* LU */

} /* prepare_vecs */


/*** update the momenta ****************************************************/
void update_anti_hermitian( site *st, int dir, Real eps, su3_matrix_f *force){
su3_matrix_f tempf1,tempf2;

    uncompress_anti_hermitian( &(st->mom[dir]), &tempf1 );
    scalar_mult_sub_su3_matrix_f( &tempf1, force, eps, &tempf2 );
    make_anti_hermitian( &tempf2, &(st->mom[dir]) );

}  /* end update_anti_hermitian */


/*** fermion force from hopping and clover terms ***************************/
void delta(field_offset psi, field_offset p)
{
  register int i, mu,nu;
  register site *st;
  msg_tag *tag0,*tag1;
  wilson_vector tvec1,tvec2;
  half_wilson_vector hvec;
  su3_matrix temp1,temp2;
  Real local_CKU0=-clov_c*kappa/(8.0*u0*u0*u0);

  for(mu=XUP;mu<=TUP;mu++){
    /**printf("mu = %d\n",mu);**/

    FORALLSITES(i,st){

      clear_su3mat(&(st->tempmat1));
      wp_shrink( (wilson_vector *)F_PT(st,psi), &(st->htmp[0]), mu, PLUS);
      wp_shrink( (wilson_vector *)F_PT(st,p), &(st->htmp[1]), mu, MINUS);
    }
    tag0 = start_gather_site( F_OFFSET(htmp[0]), sizeof(half_wilson_vector),
			      mu, EVENANDODD, gen_pt[0] );
    tag1 = start_gather_site( F_OFFSET(htmp[1]), sizeof(half_wilson_vector),
			      mu, EVENANDODD, gen_pt[1] );
    wait_gather(tag0);
    wait_gather(tag1);

    /* First the U d(dslash)/dU terms (do for only one value of nu) */
    /* SF: Fermion force due to time links at t=0 [or at t=nt-1] is zero */
    FORALLSITES(i,st) {
      /* psi and p parallel transported in from positive directions */
      mult_su3_mat_hwvec( &(st->link[mu]),
			  (half_wilson_vector *)gen_pt[0][i], &hvec);
      wp_grow( &hvec, &tvec1, mu, PLUS);
      /* i even => tvec1 = (1+gamma_mu)*U*Aodd^(-1)*D*psi,
	 i odd  => tvec1 = (1+gamma_mu)*U*psi */

      mult_su3_mat_hwvec( &(st->link[mu]),
			  (half_wilson_vector *)gen_pt[1][i], &hvec);
      wp_grow( &hvec, &tvec2, mu, MINUS);
      /* i even => tvec2 = (1-gamma_mu)*U*Aodd^(-1)*D_adj*M*psi,
	 i odd  => tvec2 = (1-gamma_mu)*U*M*psi */

      su3_projector_w( &tvec1, (wilson_vector *)F_PT(st,p), &temp1 );
      su3_projector_w( &tvec2, (wilson_vector *)F_PT(st,psi), &temp2 );
      add_su3_matrix( &temp1, &temp2, &(st->tempmat1) );
      scalar_mult_su3_matrix( &(st->tempmat1), -kappa, &(st->tempmat1) );
    } /* end loop over sites */
    cleanup_gather(tag0);
    cleanup_gather(tag1);

    /*
      SF:  the force is computed for s->link[mu].
      1. We don't update spatial links at t=0="nt", so we don't care
      what force is computed for them.
      2. In udadu_mu_nu, if mu=TUP then (s+mu)->link[nu] and (s-nu+mu)->link[nu]
      can be boundary links at t=nt.
      3. In udadu_mu_nu, if nu=TUP then (s+nu)->link[mu] can be a boundary link
      at t=nt.
    */
    /* Now the U dA/dU terms */
    for(nu=XUP;nu<=TUP;nu++) if(nu!=mu) {

      /* U dA/dU from U dM/dU */
      udadu_mu_nu( p, psi, F_OFFSET(tempmat2),
		   mu, nu, EVENANDODD );
      FORALLSITES(i,st) {
	scalar_mult_add_su3_matrix( &(st->tempmat1),
				    &(st->tempmat2), local_CKU0, &(st->tempmat1) );
      }

      /* U dA/dU from U dM^dagger/dU */
      udadu_mu_nu( psi,p, F_OFFSET(tempmat2),
		   mu, nu, EVENANDODD );
      FORALLSITES(i,st) {
	scalar_mult_add_su3_matrix( &(st->tempmat1),
				    &(st->tempmat2), local_CKU0, &(st->tempmat1) );
      }



    } /* end loop over nu & endif( nu != mu )*/
    FORALLSITES(i,st) {
      add_su3_matrix(&(st->Force[mu]),&(st->tempmat1),&(st->Force[mu]));
    }
  } /* mu */
} /* end delta */


#ifdef BETA_FREP
/*** gauge-action force from frep term **************************************/
void gauge_force_frep(field_offset scratch1, field_offset scratch2,
                      field_offset forcemat, int dir1) {
register int i,dir2;
register site *st;
msg_tag *tag0,*tag1,*tag2;
su3_matrix tmat1,tmat2;
/* divide by nflavors=2 to compensate for
ferm_epsilon = nflavors*eps = 2*eps
*/
Real bfrep_force = beta_frep/(Real)DIMF/(Real)nflavors;

/* Loop over directions for updating mom[dir1] is in the calling routine.
   Here, we loop over other directions, computing force from frep plaquettes
   in the dir1,dir2 plane
*/
	FORALLSITES(i,st){
        /* initialize staple */
            clear_su3mat( ((su3_matrix *)F_PT(st,scratch1)) );
        }
	for(dir2=XUP;dir2<=TUP;dir2++)if(dir2 != dir1){

	    /* get link[dir2] from direction dir1 */
	    tag0 = start_gather_site( F_OFFSET(link[dir2]),
                sizeof(su3_matrix), dir1, EVENANDODD, gen_pt[0] );

	    /* Start gather for the "upper staple" */
	    tag2 = start_gather_site( F_OFFSET(link[dir1]),
                sizeof(su3_matrix), dir2, EVENANDODD, gen_pt[2] );

	    /* begin the computation "at the dir2DOWN point", we will
		later gather the intermediate result "to the home point" */

	    /* Note: For SF we don't care if we get this wrong for
		     dir1<TUP and t=0, since then the staple will not be used,
		     as those links and momenta are frozen */

	    wait_gather(tag0);
            /* SF: supersede the gather with boundary links if necessary */
            FORALLSITES(i,st) {
                gen_pt[0][i] = CHOOSE_NBR(i,st,dir1,link_bndr_up[dir2],0);
            }

	    FORALLSITES(i,st){
		mult_su3_an( &(st->link[dir2]), &(st->link[dir1]), &tmat1);
		mult_su3_nn( &tmat1, (su3_matrix *)gen_pt[0][i],
		             ((su3_matrix *)F_PT(st,scratch2)) );
	    }

	    /* Gather lower staple "up to home site" */
	    tag1 = start_gather_site( scratch2, sizeof(su3_matrix),
		OPP_DIR(dir2), EVENANDODD, gen_pt[1] );

	    /*  The "upper" staple.  Note that
		one of the links has already been gathered, since it
		was used in computing the "lower" staple of the site
		above us (in dir2) */
	    wait_gather(tag2);
            /* SF: supersede the gather with boundary links if necessary */
            FORALLSITES(i,st) {
                gen_pt[2][i] = CHOOSE_NBR(i,st,dir2,link_bndr_up[dir1],2);
            }

	    FORALLSITES(i,st){
	        mult_su3_nn( &(st->link[dir2]),
       		   (su3_matrix *)gen_pt[2][i], &tmat1);
		mult_su3_na( &tmat1, (su3_matrix *)gen_pt[0][i],  &tmat2);
		add_su3_matrix( ((su3_matrix *)F_PT(st,scratch1)), &tmat2,
                                ((su3_matrix *)F_PT(st,scratch1)) );
	    }

	    wait_gather(tag1);

	    FORALLDYNLINKS(i,st,dir1){
		add_su3_matrix( ((su3_matrix *)F_PT(st,scratch1)),
                                (su3_matrix *)gen_pt[1][i],
                                ((su3_matrix *)F_PT(st,scratch1)) );
	    }
	    cleanup_gather(tag0);
	    cleanup_gather(tag1);
	    cleanup_gather(tag2);
	}

	/* Now multiply the staple sum by the link, then add to force */
	FORALLDYNLINKS(i,st,dir1){
	    mult_su3_na( &(st->link[dir1]),
                           ((su3_matrix *)F_PT(st,scratch1)), &tmat1 );
            scalar_mult_add_su3_matrix( ((su3_matrix *)F_PT(st,forcemat)),
                    &tmat1, bfrep_force, ((su3_matrix *)F_PT(st,forcemat)) );
	}
  /* no loop over dir1 here*/
} /* gauge_force_frep */
#endif /* BETA_FREP */
