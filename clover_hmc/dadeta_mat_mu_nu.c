/******************** dadeta_mat_mu_nu.c *********************************
This routine is built after udadu_mat_mu_nu with the following
main differences

1. To a force term   U dA_{oo}/dU tr_dirac(sigma_{mu,nu} A_{oo}^{-1})
   corresponds
   tr_color dU/deta dA_{oo}/dU tr_dirac(sigma_{mu,nu} A_{oo}^{-1})

2. dU/deta is non-zero only at the boundaries.
   Thus only the two lower leaves contribute at t=0,
   and only the two upper leaves contribute at t=nt.
***********************************************************************/

#include "cl_dyn_includes.h"

void dadeta_mat_mu_nu( int dir, double *f_coupling_dir ) {
register int i;
register site *s;
double f_sum;
su3_matrix mat1, mat2;
msg_tag *mtag0, *mtag1, *mtag2;

   /* get link[TUP] from direction +dir */
   mtag0 = start_gather_site( F_OFFSET(link[TUP]), sizeof(su3_matrix),
        dir, EVENANDODD, gen_pt[0] );

   /* tempmat1 = tr_dirac(gamma_4*gamma_dir*A_{oo}^{-1})  */
   tr_sigma_ldu_mu_nu_site( F_OFFSET(tempmat1), TUP, dir );

   /* get tempmat1 from direction +dir   */
   mtag1 = start_gather_site( F_OFFSET(tempmat1), sizeof(su3_matrix),
        dir, EVENANDODD, gen_pt[1] );

   FORALLSITES(i,s) {
      /* start computing staple at down boundary */
      if(s->t==0) {
	 mult_su3_an( &(s->link[TUP]), &(link_driv_dn[dir]),
               &(s->tempmat2) );
      }
   }

   wait_gather(mtag0);

   FORALLSITES(i,s) {
      if(s->t==0) {
         mult_su3_nn( &(s->tempmat2), (su3_matrix *)(gen_pt[0][i]),
               &(s->staple) );
      }
   }

   /* gather staple from TDOWN, for use at t=1 */
   mtag2 = start_gather_site( F_OFFSET(staple), sizeof(su3_matrix),
        TDOWN, EVENANDODD, gen_pt[2] );

   /* compute contribution of up boundary */
   f_sum=0.;

   /* upper-right leaf */
   FORODDSITES(i,s) {
      if(s->t==nt-1) {
         mult_su3_an( &(s->link[TUP]), &(s->tempmat1), &mat1 );
         mult_su3_nn( &mat1, &(s->link[dir]), &mat2 );
         mult_su3_nn( &mat2, (su3_matrix *)(gen_pt[0][i]), &mat1 );
         f_sum += (double)realtrace_su3( &(link_driv_up[dir]), &mat1 );
      }
   }

   wait_gather(mtag1);

   /* upper-left leaf */
   FOREVENSITES(i,s) {
      if(s->t==nt-1) {
         mult_su3_nn( &(s->link[dir]),
               (su3_matrix *)(gen_pt[1][i]), &mat1 );
         mult_su3_an( &(s->link[TUP]), &mat1, &mat2 );
         mult_su3_nn( &mat2, (su3_matrix *)(gen_pt[0][i]), &mat1 );
         f_sum += (double)realtrace_su3( &(link_driv_up[dir]), &mat1 );
      }
   }

   wait_gather(mtag2);

   /* compute contribution of down boundary */
   /* down-right leaf */
   FORODDSITES(i,s) {
      if(s->t==1) {
         mult_su3_nn( &(s->tempmat1), (su3_matrix *)(gen_pt[2][i]), &mat1 );
         f_sum += (double)realtrace_su3( &(s->link[dir]), &mat1 );
      }
   }

   /* down-left leaf */
   FOREVENSITES(i,s) {
      if(s->t==1) {
         mult_su3_nn( (su3_matrix *)(gen_pt[2][i]),
               (su3_matrix *)(gen_pt[1][i]), &mat1 );
         f_sum += (double)realtrace_su3( &(s->link[dir]), &mat1 );
      }
   }

   cleanup_gather(mtag0);
   cleanup_gather(mtag1);
   cleanup_gather(mtag2);

   g_doublesum( &f_sum );
   *f_coupling_dir = f_sum;

} /* dadeta_mat_mu_nu */
