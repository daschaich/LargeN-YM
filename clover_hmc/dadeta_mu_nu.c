/************************** coupling.c *********************************
This routine is built after udadu_mu_nu with the following
main differences

1. To a force term    U dM/dU |psi> <p|M^dag
   corresponds        <p|M^dag dU/deta dM/dU |psi>

2. dU/deta is non-zero only at the boundaries.
   Thus only the two lower leaves contribute at t=0,
   and only the two upper leaves contribute at t=nt.
***********************************************************************/

#include "cl_dyn_includes.h"

void dadeta_mu_nu( field_offset lsrc, field_offset rsrc, int dir,
                   double *f_coupling_dir ) {
register int i;
register site *s;
complex cc;
double f_sum;
wilson_vector rtemp1, rtemp2;
msg_tag *mtag0, *mtag1, *mtag2, *mtag3;

   /* get link[TUP] from direction +dir */
   mtag0 = start_gather_site( F_OFFSET(link[TUP]), sizeof(su3_matrix),
        dir, EVENANDODD, gen_pt[0] );

   FORALLSITES(i,s) {
      /* multiply rsrc by gamma_4 gamma_dir */
      if(s->t==1 || s->t==nt-1) {
         mult_sigma_mu_nu( ((wilson_vector *)F_PT(s,rsrc)),
               &(s->tmp), TUP, dir);
      }
      /* start computing staple at down boundary */
      if(s->t==0) {
	 mult_su3_an( &(s->link[TUP]), &(link_driv_dn[dir]),
               &(s->tempmat2) );
      }
   }

   /* get sources from direction +dir   */
   mtag1 = start_gather_site( lsrc, sizeof(wilson_vector),
        dir, EVENANDODD, gen_pt[1] );
   mtag2 = start_gather_site( F_OFFSET(tmp), sizeof(wilson_vector),
        dir, EVENANDODD, gen_pt[2] );

   wait_gather(mtag0);

   FORALLSITES(i,s) {
      if(s->t==0) {
         mult_su3_nn( &(s->tempmat2), (su3_matrix *)(gen_pt[0][i]),
               &(s->staple) );
      }
   }

   /* gather staple from TDOWN, for use at t=1 */
   mtag3 = start_gather_site( F_OFFSET(staple), sizeof(su3_matrix),
        TDOWN, EVENANDODD, gen_pt[3] );

   /* compute contribution of up boundary */
   f_sum=0.;

   /* upper-right leaf */
   FORALLSITES(i,s) {
      if(s->t==nt-1) {
         mult_adj_mat_wilson_vec( &(s->link[TUP]), &(s->tmp), &rtemp1 );
         mult_adj_mat_wilson_vec( &(link_driv_up[dir]),
               &rtemp1, &rtemp2 );
         mult_mat_wilson_vec( (su3_matrix *)(gen_pt[0][i]),
               &rtemp2, &rtemp1 );
         mult_mat_wilson_vec( &(s->link[dir]), &rtemp1, &rtemp2 );
         cc = wvec_dot( ((wilson_vector *)F_PT(s,lsrc)), &rtemp2 );
         f_sum += (double)cc.real;
      }
   }

   wait_gather(mtag1);
   wait_gather(mtag2);

   /* upper-left leaf */
   FORALLSITES(i,s) {
      if(s->t==nt-1) {
         mult_mat_wilson_vec( &(s->link[dir]),
               (wilson_vector *)(gen_pt[2][i]), &rtemp1 );
         mult_adj_mat_wilson_vec( &(s->link[TUP]), &rtemp1, &rtemp2 );
         mult_adj_mat_wilson_vec( &(link_driv_up[dir]),
   	       &rtemp2, &rtemp1 );
         mult_mat_wilson_vec( (su3_matrix *)(gen_pt[0][i]),
   	       &rtemp1, &rtemp2 );
         cc = wvec_dot( (wilson_vector *)(gen_pt[1][i]), &rtemp2 );
         f_sum += (double)cc.real;
      }
   }

   wait_gather(mtag3);

   /* compute contribution of down boundary */
   FORALLSITES(i,s) {
      if(s->t==1) {

         /* down-right leaf */
         mult_adj_mat_wilson_vec( &(s->link[dir]), &(s->tmp), &rtemp1 );
         mult_mat_wilson_vec( (su3_matrix *)(gen_pt[3][i]),
	       &rtemp1, &rtemp2 );
         cc = wvec_dot( ((wilson_vector *)F_PT(s,lsrc)), &rtemp2 );
         f_sum += (double)cc.real;

         /* down-left leaf */
         mult_mat_wilson_vec( (su3_matrix *)(gen_pt[3][i]),
               (wilson_vector *)(gen_pt[2][i]), &rtemp1 );
         mult_adj_mat_wilson_vec( &(s->link[dir]), &rtemp1, &rtemp2 );
         cc = wvec_dot( (wilson_vector *)(gen_pt[1][i]), &rtemp2 );
         f_sum += (double)cc.real;
      }
   }

   cleanup_gather(mtag0);
   cleanup_gather(mtag1);
   cleanup_gather(mtag2);
   cleanup_gather(mtag3);

   g_doublesum( &f_sum );
   *f_coupling_dir = f_sum;

} /* dadeta_mu_nu */
