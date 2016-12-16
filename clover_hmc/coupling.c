/************************** coupling.c *********************************
Measure the coupling in the Schroedinger functional scheme,
or more precisely K/g^2, with (for SU(3) only, and nx=ny=nz=nt,
see hep-lat/9207009 and hep-lat/9309005):

K=12*nx^2*{sin(\theta)+sin(2\theta)},   \theta=\pi/(3*nx^2)

gauge_coupling    is adapted from     generic_schroed/coupling.c
fermion_coupling  is modelled after   update_h_sf.c
dadeta_mu_nu      is modelled after   udadu_mu_nu.c
dadeta_mat_mu_nu  is modelled after   udadu_mat_mu_nu.c

printout (LU+NHYP+BETA_FREP case):

COUPL coupling bndr_plaq g_coupling bndr_plaq_frep g_coupling_frep 2*f_fat f_pseudo f_tr f_pseudo+f_tr

                                         ^            ^               ^              ^
                                     BETA_FREP     BETA_FREP         NHYP            LU

where:

bndr_plaq  = average space-time boundary plaquette

coupling = total SF coupling = K/g^2

g_coupling = -dS_plaquette/deta

f_coupling = f_pseudo [+f_tr] [+2*f_fat]
                        LU      NHYP

all cases:     f_pseudo = <psi| M^dag d(clover term)/deta |psi> + h.c.

LU:            f_tr = 2 tr (A_oo)^{-1} dA_oo/deta (again from clover term)

NHYP:          2*f_fat = <psi| M^dag dM/dV |psi> * dV/deta + h.c.
               in which V is the fat link

source: |psi> = (M^dag M)^{-1} * M^dag * g_rand

BETA_FREP:

g_coupling_frep = -dS_plaquette_frep/deta
from explicit dependence only on boundary links

Implicit dependence thru the fat links added to f_fat.

***********************************************************************/

#include "cl_dyn_includes.h"

void gauge_coupling(double *bndr_plaq, double *g_coupling);
#ifdef BETA_FREP
void gauge_coupling_frep(double *bndr_plaq_frep, double *g_coupling_frep);
#endif
#ifdef LU
void fermion_coupling(double *f_pseudo, double *f_tr);
#else
void fermion_coupling(double *f_pseudo);
#endif

/*** calculate SF's K/g^2 *********************************************/
void coupling() {
double bndr_plaq, g_coupling, f_pseudo;
#ifdef BETA_FREP
double bndr_plaq_frep, g_coupling_frep;
#endif
#ifdef LU
double f_tr;
#endif

    /* gauge field part of 1/g^2              */
    gauge_coupling(&bndr_plaq, &g_coupling);
#ifdef BETA_FREP
    gauge_coupling_frep(&bndr_plaq_frep, &g_coupling_frep);
#endif
    /* fermion part of 1/g^2 and printout     */
    /* total K/g^2 is always the first output */

#ifdef LU

    fermion_coupling(&f_pseudo, &f_tr);

#ifdef NHYP
#ifdef BETA_FREP
    if(this_node==0)printf("COUPL %e %e %e %e %e %e %e %e %e\n",
        g_coupling+g_coupling_frep+f_pseudo+f_tr+2*f_fat,
        bndr_plaq, g_coupling,
        bndr_plaq_frep, g_coupling_frep,
        2*f_fat, f_pseudo, f_tr, f_pseudo+f_tr );
#else
    if(this_node==0)printf("COUPL %e %e %e %e %e %e %e\n",
        g_coupling+f_pseudo+f_tr+2*f_fat,
        bndr_plaq, g_coupling,
        2*f_fat, f_pseudo, f_tr, f_pseudo+f_tr );
#endif /* BETA_FREP */
#else  /* thin links */
#ifdef BETA_FREP
    if(this_node==0)printf("COUPL %e %e %e %e %e %e %e %e\n",
        g_coupling+g_coupling_frep+f_pseudo+f_tr,
        bndr_plaq, g_coupling,
        bndr_plaq_frep, g_coupling_frep,
        f_pseudo, f_tr, f_pseudo+f_tr );
#else
    if(this_node==0)printf("COUPL %e %e %e %e %e %e\n",
        g_coupling+f_pseudo+f_tr,
        bndr_plaq, g_coupling,
        f_pseudo, f_tr, f_pseudo+f_tr );
#endif /* BETA_FREP */
#endif /* NHYP  */

#else  /* no LU */

    fermion_coupling(&f_pseudo);

#ifdef NHYP
#ifdef BETA_FREP
    if(this_node==0)printf("COUPL %e %e %e %e %e %e %e\n",
        g_coupling+g_coupling_frep+f_pseudo+2*f_fat,
        bndr_plaq, g_coupling,
        bndr_plaq_frep, g_coupling_frep,
        2*f_fat, f_pseudo );
#else
    if(this_node==0)printf("COUPL %e %e %e %e %e\n",
        g_coupling+f_pseudo+2*f_fat,
        bndr_plaq, g_coupling,
        2*f_fat, f_pseudo );
#endif /* BETA_FREP */
#else  /* thin links */
#ifdef BETA_FREP
    if(this_node==0)printf("COUPL %e %e %e %e %e %e\n",
        g_coupling+g_coupling_frep+f_pseudo,
        bndr_plaq, g_coupling,
        bndr_plaq_frep, g_coupling_frep,
        f_pseudo );
#else
    if(this_node==0)printf("COUPL %e %e %e %e\n",
        g_coupling+f_pseudo,
        bndr_plaq, g_coupling,
        f_pseudo );
#endif /* BETA_FREP */
#endif /* NHYP  */

#endif /* LU    */

} /* coupling */


/***  Gauge action contribution to 1/g^2 ******************************/
void gauge_coupling(double *bndr_plaq, double *g_coupling) {
register int i,dir;
register site *s;
su3_matrix_f mtmp;
double bndr_sum=0., driv_sum=0.;
msg_tag *mtag0, *mtag1;

   for(dir=XUP;dir<=ZUP;dir++){

      mtag0 = start_gather_site( F_OFFSET(linkf[dir]), sizeof(su3_matrix_f),
	    TUP, EVENANDODD, gen_pt[0] );
      mtag1 = start_gather_site( F_OFFSET(linkf[TUP]), sizeof(su3_matrix_f),
	    dir, EVENANDODD, gen_pt[1] );

      FORALLSITES(i,s){
         if(s->t == 0){
	    mult_su3_an_f( &(s->linkf[TUP]), &(linkf_bndr_dn[dir]),
               (su3_matrix_f *)&(s->tempmat1) );
	    mult_su3_an_f( &(s->linkf[TUP]), &(linkf_driv_dn[dir]),
               (su3_matrix_f *)&(s->tempmat2) );
         }
	 else if(s->t == nt-1){
	    mult_su3_an_f( &(s->linkf[TUP]), &(s->linkf[dir]),
               (su3_matrix_f *)&(s->staple) );
         }
      }

      wait_gather(mtag0);
      wait_gather(mtag1);

      FORALLSITES(i,s){
         if(s->t == 0){
	    mult_su3_nn_f( (su3_matrix_f *)&(s->tempmat1),
	       (su3_matrix_f *)(gen_pt[1][i]), &mtmp);
	    bndr_sum += (double)
	       realtrace_su3_f( (su3_matrix_f *)(gen_pt[0][i]), &mtmp );
	    mult_su3_nn_f( (su3_matrix_f *)&(s->tempmat2),
	       (su3_matrix_f *)(gen_pt[1][i]), &mtmp);
	    driv_sum += (double)
	       realtrace_su3_f( (su3_matrix_f *)(gen_pt[0][i]), &mtmp );
         }
	 else if(s->t == nt-1){
	    mult_su3_nn_f( (su3_matrix_f *)&(s->staple),
	       (su3_matrix_f *)(gen_pt[1][i]), &mtmp);
	    bndr_sum += (double)
	       realtrace_su3_f( &(linkf_bndr_up[dir]), &mtmp );
	    driv_sum += (double)
	       realtrace_su3_f( &(linkf_driv_up[dir]), &mtmp );
         }
      }

   cleanup_gather(mtag0);
   cleanup_gather(mtag1);
   } /* loop over dir */

   g_doublesum( &bndr_sum );
   g_doublesum( &driv_sum );

   *bndr_plaq = bndr_sum/((double)(6*nx*ny*nz));
   *g_coupling = -driv_sum*((double)beta)/((double)NCOL);
} /* gauge_coupling */


/***  BETA_FREP Gauge action contribution to 1/g^2 ********************/
#ifdef BETA_FREP
void gauge_coupling_frep(double *bndr_plaq_frep, double *g_coupling_frep) {
register int i,dir;
register site *s;
su3_matrix mtmp;
double bndr_sum_frep=0., driv_sum_frep=0.;
msg_tag *mtag0, *mtag1;

   for(dir=XUP;dir<=ZUP;dir++){

      mtag0 = start_gather_site( F_OFFSET(link[dir]), sizeof(su3_matrix),
	    TUP, EVENANDODD, gen_pt[0] );
      mtag1 = start_gather_site( F_OFFSET(link[TUP]), sizeof(su3_matrix),
	    dir, EVENANDODD, gen_pt[1] );

      FORALLSITES(i,s){
         if(s->t == 0){
	    mult_su3_an( &(s->link[TUP]), &(link_bndr_dn[dir]),
                         &(s->tempmat1) );
	    mult_su3_an( &(s->link[TUP]), &(link_driv_dn[dir]),
                         &(s->tempmat2) );
         }
	 else if(s->t == nt-1){
	    mult_su3_an( &(s->link[TUP]), &(s->link[dir]),&(s->staple) );
         }
      }

      wait_gather(mtag0);
      wait_gather(mtag1);

      FORALLSITES(i,s){
         if(s->t == 0){
	    mult_su3_nn( &(s->tempmat1), (su3_matrix *)(gen_pt[1][i]), &mtmp);
	    bndr_sum_frep += (double)
	       realtrace_su3( (su3_matrix *)(gen_pt[0][i]), &mtmp );
	    mult_su3_nn( &(s->tempmat2), (su3_matrix *)(gen_pt[1][i]), &mtmp);
	    driv_sum_frep += (double)
	       realtrace_su3( (su3_matrix *)(gen_pt[0][i]), &mtmp );
         }
	 else if(s->t == nt-1){
	    mult_su3_nn( &(s->staple), (su3_matrix *)(gen_pt[1][i]), &mtmp);
	    bndr_sum_frep += (double)
	       realtrace_su3( &(link_bndr_up[dir]), &mtmp );
	    driv_sum_frep += (double)
	       realtrace_su3( &(link_driv_up[dir]), &mtmp );
         }
      }

   cleanup_gather(mtag0);
   cleanup_gather(mtag1);
   } /* loop over dir */

   g_doublesum( &bndr_sum_frep );
   g_doublesum( &driv_sum_frep );

   *bndr_plaq_frep = bndr_sum_frep/((double)(6*nx*ny*nz));
   *g_coupling_frep = -driv_sum_frep*((double)beta_frep)/((double)DIMF);
} /* gauge_coupling_frep */
#endif /* BETA_FREP */


/*** fermion's contribution to 1/g^2 **********************************/
#ifdef LU
void fermion_coupling(double *f_pseudo, double *f_tr) {
#else
void fermion_coupling(double *f_pseudo) {
#endif

register int i, dir;
register site *s;

int num_masses_tmp;

int iters;
Real final_rsq;
double f_coupling_dir;
Real CKU0=clov_c*kappa/(u0*u0*u0);
  iters=0;
   make_clov(CKU0);
#ifdef LU
   make_clovinv(ODD);
#endif

  /* kluge to avoid hasenbusch */
	num_masses_tmp=num_masses;
	num_masses=1;


   /* get a new random source: if num_masses=1, this will appear from grsource_w() in chi[0].
      We also need the fermion force. Again, only if num_masses=1 will this be computed
      using the unshifted, unreweighted valence quark mass
   chi = M^dag * g_rand
   psi = (M^dag M)^{-1} * chi -- note inversion does NOT use the ``shift'' parameter!
   LU: all live on even sites only  */
   grsource_w();

   /* zero initial guess */
   FORALLSITESDOMAIN(i,s) {
      clear_wvec( &(s->psi[0]) );
   }

  iters += congrad_cl_m(niter,rsqmin,&final_rsq,F_OFFSET(chi[0]),F_OFFSET(psi[0]),0.0);

#ifdef NHYP
   /* fat links.  fermion_force will do the inversion,
      compute the fat-links Sigma, and call stout_force1
      to compute the contribution from the dependence of
      the fat links on the boundary links                          */

   /* initialize globals before calling fermion_force              */
   /* no momentum update anyway, so we set eps to arbitrary value  */
   sf_coupling_flag=SF_COUPLING;
   f_fat=0.;
   fermion_force(0.0,0.0);
   g_doublesum( &f_fat );
#else
   /* thin links.  we only need to do the inversion                */
   prepare_vecs(0);
#endif

   /* now we can restore the program's parameters */
  num_masses=num_masses_tmp;


   /* now the clovers' contribution to K/g^2                       */
   *f_pseudo = 0.;
#ifdef LU
   *f_tr = 0.;
#endif

   for(dir=XUP;dir<=ZUP;dir++){
      /* add: dA/deta from <p|dM/deta|psi>     */
      dadeta_mu_nu( F_OFFSET(p), F_OFFSET(psi[0]), dir, &f_coupling_dir );
      *f_pseudo += f_coupling_dir;
      /* add: dA/deta from <psi|dM^dag/deta|p> */
      dadeta_mu_nu( F_OFFSET(psi[0]), F_OFFSET(p), dir, &f_coupling_dir );
      *f_pseudo += f_coupling_dir;

#ifdef LU
      /* add:
      tr_color dU/deta dA_{oo}/dU tr_dirac(sigma_{mu,nu} A_{oo}^{-1})
      */
      dadeta_mat_mu_nu( dir, &f_coupling_dir );
      *f_tr += f_coupling_dir;
#endif /*LU*/

   } /* end loop over dir */

   /* normalization and sign*/
   *f_pseudo *= (double)(-CKU0/4.);
#ifdef LU
   *f_tr *= (double)(-CKU0/2.);
#endif

   free_clov();

} /* fermion_coupling */


