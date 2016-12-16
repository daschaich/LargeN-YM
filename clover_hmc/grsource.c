/* Clover fermions */

#include "cl_dyn_includes.h"

/* two-mass version. if num_masses=1, then usual code. If 2, then
chi[0] and psi[0] are the ``slow'' (no-shift) inverter and chi[1] and psi[1]
are the fast (shifted) pseudofermions */

/* construct a gaussian random vector, g_rand, and chi=M(dagger)*g_rand  */
/* also clear psi, since zero is our best guess for the solution with a
   new random chi field. */

/* if "LU" is defined use the LU preconditioned fermion matrix, where
  the fermion spinors live on even sites only.  In other words, if
  Dslash_oe is the dslash operator with its source on even sites and
  its result on odd sites, etc.:

  without LU:
   M = A - kappa*( Dslash_eo + Dslash_oe )
  with LU:
   M = A_even - kappa^2 * Dslash_eo * (A_odd)^{-1} * Dslash_oe
*/
#ifdef LU
#define FORMYSITESDOMAIN FOREVENSITESDOMAIN
#else
#define FORMYSITESDOMAIN FORALLSITESDOMAIN
#endif

void grsource_w() {
register int i,j,k;
register site *s;
 int kind1;
 wilson_vector tmpvec,tmpvec2;
 double rr;
 int iters;
 Real final_rsq;

  iters=0;

/* First the original field or the simply shifted one in the 2PF case */
 if(num_masses==1) kind1=0;
 else kind1=1;


 rr=0.0;
 FORMYSITESDOMAIN(i,s){
   for(k=0;k<4;k++)for(j=0;j<DIMF;j++){
#ifdef SITERAND
     s->g_rand.d[k].c[j].real = gaussian_rand_no(&(s->site_prn));
     s->g_rand.d[k].c[j].imag = gaussian_rand_no(&(s->site_prn));
#else
     s->g_rand.d[k].c[j].real = gaussian_rand_no(&node_prn);
     s->g_rand.d[k].c[j].imag = gaussian_rand_no(&node_prn);
#endif
     s->psi[kind1].d[k].c[j] = cmplx(0.0,0.0);
   }
   rr += (double)wvec_rdot(&(s->g_rand),&(s->g_rand));
 }

 
#ifdef LU
 /*  chi <- M^dagger g_rand  */
 mult_ldu_site( F_OFFSET(g_rand), F_OFFSET(tmp), EVEN );
 dslash_w_site( F_OFFSET(g_rand), F_OFFSET(tmp), MINUS, ODD);
 mult_ldu_site( F_OFFSET(tmp), F_OFFSET(g_rand), ODD );
 dslash_w_site( F_OFFSET(g_rand), F_OFFSET(chi[kind1]), MINUS, EVEN);
 FOREVENSITESDOMAIN(i,s){
   if(num_masses==2){
     scalar_mult_add_wvec( &(s->tmp), &(s->chi[kind1]), -kappa*kappa,
			   &tmpvec );
     /* That was Mdagger, now we need to subtract -i*shift*gamma5*p */
     mult_by_gamma( &(s->g_rand), &tmpvec2, GAMMAFIVE );
     c_scalar_mult_add_wvec2( &tmpvec, &tmpvec2, cmplx(0.0,-shift), &(s->chi[kind1]) );
   }
   else{
     scalar_mult_add_wvec( &(s->tmp), &(s->chi[kind1]), -kappa*kappa,
			   &(s->chi[kind1]) );
   }
 }
#else
 mult_ldu_site( F_OFFSET(g_rand), F_OFFSET(tmp), EVENANDODD );
 dslash_w_site( F_OFFSET(g_rand), F_OFFSET(chi[kind1]), MINUS, EVENANDODD);
 FORALLSITESDOMAIN(i,s){
   if(num_masses==2){
     scalar_mult_add_wvec( &(s->tmp), &(s->chi[kind1]), -kappa,
			   &(tmpvec) );
     mult_by_gamma( &(s->g_rand), &tmpvec2, GAMMAFIVE );
     c_scalar_mult_add_wvec2( &tmpvec, &tmpvec2, cmplx(0.0,-shift), &(s->chi[kind1]) );
   }
     else{
       scalar_mult_add_wvec( &(s->tmp), &(s->chi[kind1]), -kappa,
			     &(s->chi[kind1]) );
     }
 }
#endif
 
 /* Now the slightly more fancy case of the 2nd PF field: chi[0] and psi[0] */
 if(num_masses==2){
   /*
   node0_printf("CALL 2nd PF\n");
   */
   /* The first contribution is just the random number itself
      so we might as well generate it in chi[0] directly */
   FORMYSITESDOMAIN(i,s){
     for(k=0;k<4;k++)for(j=0;j<DIMF;j++){
       /* Comment these lines out, set s->chi[0].d[k].c[j] = cmplx(1.0,0.0);
	  and set ``step=0.0'' in the input file and the 1 and 2 PF codes should produce
	  identical results */
#ifdef SITERAND
       s->chi[0].d[k].c[j].real = gaussian_rand_no(&(s->site_prn));
       s->chi[0].d[k].c[j].imag = gaussian_rand_no(&(s->site_prn));
#else
       s->chi[0].d[k].c[j].real = gaussian_rand_no(&node_prn);
       s->chi[0].d[k].c[j].imag = gaussian_rand_no(&node_prn);
#endif

       s->psi[0].d[k].c[j] = cmplx(0.0,0.0);
     }
     rr += (double)wvec_rdot(&(s->chi[0]),&(s->chi[0]));
   }
   
   /* But now something fancy has to happen. First we invert with the shift
      and write the result into old_psi[0] */
   iters+=congrad_cl_m(niter,rsqmin,&final_rsq, F_OFFSET(chi[0]), F_OFFSET(old_psi[0]), shift);
   
#ifdef LU
   /*  chi <- Mtilde old_psi1  */
   /*  chi <- M^dagger chi[0]  */
   mult_ldu_site( F_OFFSET(old_psi[0]), F_OFFSET(tmp), EVEN );
   dslash_w_site( F_OFFSET(old_psi[0]), F_OFFSET(tmp), PLUS, ODD);
   mult_ldu_site( F_OFFSET(tmp), F_OFFSET(old_psi[0]), ODD );
   dslash_w_site( F_OFFSET(old_psi[0]), F_OFFSET(p), PLUS, EVEN);
   FOREVENSITESDOMAIN(i,s){
     scalar_mult_add_wvec( &(s->tmp), &(s->p), -kappa*kappa, &tmpvec);
     
     /* That was M, now we need to add i*shift*gamma5*p */
     mult_by_gamma( &(s->old_psi[0]), &tmpvec2, GAMMAFIVE );
     c_scalar_mult_add_wvec2( &tmpvec, &tmpvec2, cmplx(0.0,shift), &tmpvec );
     /* And the final step is to multiply by i*shift*gamma5*
	before adding it to chi[0] */
     mult_by_gamma( &tmpvec, &tmpvec2, GAMMAFIVE );
     c_scalar_mult_add_wvec2( &(s->chi[0]), &tmpvec2, cmplx(0.0,shift), &(s->chi[0]));
   }
#else
   mult_ldu_site( F_OFFSET(old_psi[0]), F_OFFSET(tmp), EVENANDODD );
   dslash_w_site( F_OFFSET(old_psi[0]), F_OFFSET(p), PLUS, EVENANDODD);
   FORALLSITESDOMAIN(i,s){
     scalar_mult_add_wvec( &(s->tmp), &(s->p), -kappa,&(tmpvec) );
     /* That was M, now we need to add i*shift*gamma5*old_psi[0] */
     mult_by_gamma( &(s->old_psi[0]), &tmpvec2, GAMMAFIVE );
     c_scalar_mult_add_wvec2( &tmpvec, &tmpvec2, cmplx(0.0,shift), &tmpvec );
     /* And the final step is to multiply by i*shift*gamma5*
	before adding it to chi2 */
     mult_by_gamma( &tmpvec, &tmpvec2, GAMMAFIVE );
     c_scalar_mult_add_wvec2( &(s->chi[0]), &tmpvec2, cmplx(0.0,shift), &(s->chi[0]));
   }
   
   
   
#endif
 
 }
 
 g_doublesum(&rr);
 /*
 node0_printf("rr %e\n",rr);
 */  
 
}/* grsource_w */


/* Check congrad by multiplying psi by Madj*M, compare result to chi */
/* Before calling checkmul() you should call  congrad_w() */
void checkmul(field_offset chi, field_offset psi, Real mshift) {
register int i,j,k;
wilson_vector wvec;
register site *s;

 printf("CHECKMUL: starting %e\n",mshift);

    /* multiply by M_adjoint*M */
#ifdef LU
    mult_ldu_site( psi, F_OFFSET(tmp), EVEN );
    dslash_w_site( psi, psi, PLUS, ODD );
    mult_ldu_site( psi, F_OFFSET(mp),ODD );
    dslash_w_site( F_OFFSET(mp), F_OFFSET(mp), PLUS, EVEN);
    FOREVENSITESDOMAIN(i,s)
        scalar_mult_add_wvec( &(s->tmp), &(s->mp), -kappa*kappa, &(s->mp) );
    mult_ldu_site(F_OFFSET(mp), F_OFFSET(tmp), EVEN);
    dslash_w_site( F_OFFSET(mp), F_OFFSET(mp), MINUS, ODD );
    mult_ldu_site(F_OFFSET(mp), F_OFFSET(tmp), ODD);
    dslash_w_site( F_OFFSET(tmp), F_OFFSET(p) , MINUS, EVEN);
    FOREVENSITESDOMAIN(i,s){
        scalar_mult_add_wvec( &(s->tmp), &(s->p), -kappa*kappa, &(s->p));
	scalar_mult_add_wvec(&(s->p),(wilson_vector *)F_PT(s,psi),mshift*mshift,&(s->p));
   }
#else
    mult_ldu_site( psi), F_OFFSET(tmp), EVENANDODD );
    dslash_w_site( psi, F_OFFSET(mp), PLUS, EVENANDODD);
    FORALLSITESDOMAIN(i,s)
        scalar_mult_add_wvec( &(s->tmp), &(s->mp), -kappa, &(s->mp) );
    mult_ldu_site( F_OFFSET(mp), F_OFFSET(tmp), EVENANDODD );
    dslash_w_site( F_OFFSET(mp), F_OFFSET(p), MINUS, EVENANDODD);
    FORALLSITESDOMAIN(i,s){
        scalar_mult_add_wvec( &(s->tmp), &(s->p), -kappa, &(s->p) );
       scalar_mult_add_wvec(&(s->p),(wilson_vector *)F_PT(s,psi),mshift*mshift,&(s->p));
   }
#endif

    /* Compare to source */
    FORMYSITESDOMAIN(i,s){
	/**printf("Site %d %d %d %d\n",s->x,s->y,s->z,s->t);**/
      copy_wvec((wilson_vector *)F_PT(s,chi),&wvec);
	for(k=0;k<4;k++)for(j=0;j<DIMF;j++){

	  if( fabs((double)(wvec.d[k].c[j].real  -(double)s->p.d[k].c[j].real) > 1e-4 )){
	    printf("%d %d %d\t%e\t%e\t%e\n",i,k,j,
(double)(wvec.d[k].c[j].real),
		(double)s->p.d[k].c[j].real,
(double)(wvec.d[k].c[j].real) - (double)s->p.d[k].c[j].real);
if( fabs((double)(wvec.d[k].c[j].imag)-(double)s->p.d[k].c[j].imag) > 1e-4 )
	    printf("%d %d %d\t%e\t%e\t%e\n",i,k,j,
(double)(wvec.d[k].c[j].imag),
		(double)s->p.d[k].c[j].imag,
(double)(wvec.d[k].c[j].imag) -  (double)s->p.d[k].c[j].imag);
	}
	  /*
	printf("\n");
	  */
	}
    }
}/* checkmul */
