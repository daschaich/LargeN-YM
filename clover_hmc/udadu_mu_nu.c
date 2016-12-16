/******************* udadu_mu_nu.c ****************************************/

/* MIMD version 6 */
/* This version uses gathers to get the neighbors */
/* version of April 1995 by MBW */

/* Compute U dA/dU, part of the fermion force term, given mu & nu.

   The udadul terms are left-multiplied by the transpose of a
   wilson_vector, so we take the hermitian conjugate here, so we
   can operate on the vector from the left.  Essentially:
     x^* A = (A^\dagger x)^*
   The udadur terms are right-multiplied by a wilson_vector.

   The force term consists of only the right half of the "clover".

   lsrc and rsrc are input wilson_vectors.
   mat is output: the contribution to iH_dot
   Uses tmp as an internal wilson_vector.  psi, chi, p and mp are busy.
*/

#include "cl_dyn_includes.h"

void udadu_mu_nu( field_offset lsrc, field_offset rsrc, field_offset mat,
		 int mu, int nu, int parity )
{
register int i, otherparity = 0;
int disp[4];
register site *s;
msg_tag *tag[8];
wilson_vector ltemp, rtemp;		/* input to su3_projector_w() */
wilson_vector sittemp;
/* SF: sittemp avoids using s->tmp, which must be kept zero at t=0
 for correct functionality. */
su3_matrix tmpmat;			/* for middle steps */

    switch(parity) {
        case EVEN:              otherparity = ODD; break;
        case ODD:               otherparity = EVEN; break;
        case EVENANDODD:        otherparity = EVENANDODD; break;
    }

/**********  First work on upper leaf  **********/
    /* get link[nu] from direction +mu */
    tag[0] = start_gather_site( F_OFFSET(link[nu]), sizeof(su3_matrix),
        mu, EVENANDODD, gen_pt[0] );

    /* get link[mu] from direction +nu */
    tag[1] = start_gather_site( F_OFFSET(link[mu]), sizeof(su3_matrix),
        nu, EVENANDODD, gen_pt[1] );

 /*****  Upper leaf - same parity *****/
    /*get sources from +mu +nu */
    for(i=XUP;i<=TUP;i++) disp[i]=0;
    disp[mu] = 1;  disp[nu] = 1;
    tag[2] = start_general_gather_site( lsrc, sizeof(wilson_vector),
        disp, parity, gen_pt[2] );

    wait_gather(tag[0]);
    wait_gather(tag[1]);
/* SF: supersede the gathers with boundary links if necessary */
    FORALLSITES(i,s) {
        gen_pt[0][i] = CHOOSE_NBR(i,s,mu,link_bndr_up[nu],0);
        gen_pt[1][i] = CHOOSE_NBR(i,s,nu,link_bndr_up[mu],1);
    }
    wait_general_gather(tag[2]);
    tag[3] = start_general_gather_site( rsrc, sizeof(wilson_vector),
        disp, parity, gen_pt[3] );

    /* ltemp = lsrc, rtemp = plaq*rsrc */
    FORSOMEPARITY(i,s,parity) {
	mult_adj_mat_wilson_vec( &(s->link[nu]),
	    ((wilson_vector *)F_PT(s,rsrc)), &rtemp );
	mult_adj_mat_wilson_vec( (su3_matrix *)(gen_pt[1][i]), &rtemp,
	    &sittemp );
	mult_mat_wilson_vec( (su3_matrix *)(gen_pt[0][i]),
	    &sittemp, &rtemp );
	mult_mat_wilson_vec( &(s->link[mu]), &rtemp, &sittemp );
	mult_sigma_mu_nu( &sittemp, &rtemp, mu, nu );

	su3_projector_w( &rtemp, ((wilson_vector *)F_PT(s,lsrc)),
	    ((su3_matrix *)F_PT(s,mat)) );
    }

    wait_general_gather(tag[3]);

    /* ltemp = upperleftcorner*lsrc, rtemp = lowerrightcorner*rsrc */
    FORSOMEPARITY(i,s,parity) {
	mult_mat_wilson_vec( (su3_matrix *)(gen_pt[1][i]),
	    (wilson_vector *)(gen_pt[2][i]), &sittemp );
	mult_mat_wilson_vec( &(s->link[nu]), &sittemp, &ltemp );

	mult_mat_wilson_vec( (su3_matrix *)(gen_pt[0][i]),
	    (wilson_vector *)(gen_pt[3][i]), &rtemp );
	mult_mat_wilson_vec( &(s->link[mu]), &rtemp, &sittemp );
	mult_sigma_mu_nu( &sittemp, &rtemp, mu, nu );

	su3_projector_w( &rtemp, &ltemp, &tmpmat );
	add_su3_matrix( ((su3_matrix *)F_PT(s,mat)),
	    &tmpmat, ((su3_matrix *)F_PT(s,mat)) );
    }

 /***** upper leaf - different parity *****/

    /* get sources from +nu */
    tag[4] = start_gather_site( lsrc, sizeof(wilson_vector),
        nu, otherparity, gen_pt[4] );
    tag[5] = start_gather_site( rsrc, sizeof(wilson_vector),
        nu, otherparity, gen_pt[5] );

    wait_gather(tag[4]);
    wait_gather(tag[5]);

    /* ltemp = link[nu]*lsrc, rtemp = staple*rsrc */
    FORSOMEPARITY(i,s,otherparity) {
	mult_mat_wilson_vec( &(s->link[nu]),
	    (wilson_vector *)(gen_pt[4][i]), &ltemp );

	mult_adj_mat_wilson_vec( (su3_matrix *)(gen_pt[1][i]),
	    (wilson_vector *)(gen_pt[5][i]), &sittemp );
	mult_mat_wilson_vec( (su3_matrix *)(gen_pt[0][i]), &sittemp, &rtemp);
	mult_mat_wilson_vec( &(s->link[mu]), &rtemp, &sittemp );
	mult_sigma_mu_nu( &sittemp, &rtemp, mu, nu );

        if(parity==EVENANDODD) {
            su3_projector_w( &rtemp, &ltemp, &tmpmat );
            add_su3_matrix( ((su3_matrix *)F_PT(s,mat)),
	        &tmpmat, ((su3_matrix *)F_PT(s,mat)) );
        }
        else {
    /* mat hasn't been started on opposite parity sites yet */
            su3_projector_w( &rtemp, &ltemp, ((su3_matrix *)F_PT(s,mat)) );
        }
    }

    cleanup_gather(tag[4]);
    cleanup_gather(tag[5]);

    /* get sources from +mu */
    tag[4] = start_gather_site( lsrc, sizeof(wilson_vector),
        mu, otherparity, gen_pt[4] );
    tag[5] = start_gather_site( rsrc, sizeof(wilson_vector),
        mu, otherparity, gen_pt[5] );

    wait_gather(tag[4]);
    wait_gather(tag[5]);

    /* ltemp = staple*lsrc, rtemp = link[mu]*rsrc */
    FORSOMEPARITY(i,s,otherparity) {
	mult_adj_mat_wilson_vec( (su3_matrix *)(gen_pt[0][i]),
	    (wilson_vector *)(gen_pt[4][i]), &ltemp );
	mult_mat_wilson_vec( (su3_matrix *)(gen_pt[1][i]), &ltemp,
	    &sittemp );
	mult_mat_wilson_vec( &(s->link[nu]), &sittemp, &ltemp );

	mult_mat_wilson_vec( &(s->link[mu]),
	    (wilson_vector *)(gen_pt[5][i]), &sittemp );
	mult_sigma_mu_nu( &sittemp, &rtemp, mu, nu );

	su3_projector_w( &rtemp, &ltemp, &tmpmat );
	add_su3_matrix( ((su3_matrix *)F_PT(s,mat)),
	    &tmpmat, ((su3_matrix *)F_PT(s,mat)) );
    }

    /* Hold onto these these sources since they appear in lower leaf below */

/**********  The Lower leaf  **********/

    cleanup_gather(tag[0]);
    cleanup_gather(tag[1]);
    /* get link[nu] from direction -nu */
    tag[0] = start_gather_site( F_OFFSET(link[nu]), sizeof(su3_matrix),
        OPP_DIR(nu), EVENANDODD, gen_pt[0] );

    /* get link[mu] from direction -nu */
    tag[1] = start_gather_site( F_OFFSET(link[mu]), sizeof(su3_matrix),
        OPP_DIR(nu), EVENANDODD, gen_pt[1] );

    /* get link[nu] from direction -nu +mu */
    /* disp[mu] = 1; already */
    disp[nu] = -1;
    tag[6] = start_general_gather_site( F_OFFSET(link[nu]), sizeof(su3_matrix),
        disp, EVENANDODD, gen_pt[6] );

    wait_gather(tag[0]);
    wait_gather(tag[1]);
    wait_general_gather(tag[6]);

/* SF: supersede the general gather with boundary links if necessary.
The assignment of boundary link is (s+mu)->link[nu] =  (s+mu-nu)->link[nu],
using that boundary links are spatially constant.
(Note that (s-nu)->link[nu] and (s-nu)->link[mu] are ok.)
*/
    FORALLSITES(i,s) {
        gen_pt[6][i] = CHOOSE_NBR(i,s,mu,link_bndr_up[nu],6);
    }

 /***** lower leaf, different parity *****/

    /* ltemp = link[mu]*lsrc, rtemp = staple*rsrc */
    FORSOMEPARITY(i,s,otherparity) {
	mult_mat_wilson_vec( &(s->link[mu]),
	    (wilson_vector *)(gen_pt[4][i]), &ltemp );

	mult_mat_wilson_vec( (su3_matrix *)(gen_pt[6][i]),
	    (wilson_vector *)(gen_pt[5][i]), &sittemp );
	mult_mat_wilson_vec( (su3_matrix *)(gen_pt[1][i]),
	    &sittemp, &rtemp );
	mult_adj_mat_wilson_vec( (su3_matrix *)(gen_pt[0][i]),
	    &rtemp, &sittemp );
	mult_sigma_mu_nu( &sittemp, &rtemp, mu, nu );

	su3_projector_w( &rtemp, &ltemp, &tmpmat );
	sub_su3_matrix( ((su3_matrix *)F_PT(s,mat)),
	    &tmpmat, ((su3_matrix *)F_PT(s,mat)) );
    }

    cleanup_gather(tag[4]);
    cleanup_gather(tag[5]);
    /* get sources from -nu */
    tag[4] = start_gather_site( lsrc, sizeof(wilson_vector),
        OPP_DIR(nu), otherparity, gen_pt[4] );
    tag[5] = start_gather_site( rsrc, sizeof(wilson_vector),
        OPP_DIR(nu), otherparity, gen_pt[5] );

    wait_gather(tag[4]);
    wait_gather(tag[5]);

    /* ltemp = staple*lsrc, rtemp = link[nu]^dagger*rsrc */
    FORSOMEPARITY(i,s,otherparity) {
	mult_adj_mat_wilson_vec( (su3_matrix *)(gen_pt[1][i]),
	    (wilson_vector *)(gen_pt[4][i]), &ltemp );
	mult_adj_mat_wilson_vec( (su3_matrix *)(gen_pt[6][i]),
	    &ltemp, &sittemp );
	mult_mat_wilson_vec( &(s->link[mu]), &sittemp, &ltemp );

	mult_adj_mat_wilson_vec( (su3_matrix *)(gen_pt[0][i]),
	    (wilson_vector *)(gen_pt[5][i]), &sittemp );
	mult_sigma_mu_nu( &sittemp, &rtemp, mu, nu );

	su3_projector_w( &rtemp, &ltemp, &tmpmat );
	sub_su3_matrix( ((su3_matrix *)F_PT(s,mat)),
	    &tmpmat, ((su3_matrix *)F_PT(s,mat)) );
    }

 /***** lower leaf, same parity *****/

    cleanup_general_gather(tag[2]);
    cleanup_general_gather(tag[3]);
    /* get sources from +mu -nu */
    /* disp[mu] = 1;  disp[nu] = -1; already */
    tag[2] = start_general_gather_site( lsrc, sizeof(wilson_vector),
        disp, parity, gen_pt[2] );
    wait_general_gather(tag[2]);
    tag[3] = start_general_gather_site( rsrc, sizeof(wilson_vector),
        disp, parity, gen_pt[3] );

    /* ltemp = plaq*lsrc, rtemp = rsrc */
    FORSOMEPARITY(i,s,parity) {
	mult_mat_wilson_vec( (su3_matrix *)(gen_pt[0][i]),
	    ((wilson_vector *)F_PT(s,lsrc)), &sittemp );
	mult_adj_mat_wilson_vec( (su3_matrix *)(gen_pt[1][i]),
	    &sittemp, &ltemp );
	mult_adj_mat_wilson_vec( (su3_matrix *)(gen_pt[6][i]),
	    &ltemp, &sittemp );
	mult_mat_wilson_vec( &(s->link[mu]), &sittemp, &ltemp );

	mult_sigma_mu_nu( ((wilson_vector *)F_PT(s,rsrc)), &rtemp, mu, nu );

	su3_projector_w( &rtemp, &ltemp, &tmpmat );
	sub_su3_matrix( ((su3_matrix *)F_PT(s,mat)),
	    &tmpmat, ((su3_matrix *)F_PT(s,mat)) );
    }

    wait_general_gather(tag[3]);

    /* ltemp = upperrightcorner*lsrc, rtemp = lowerleftcorner*rsrc */
    FORSOMEPARITY(i,s,parity) {
	mult_adj_mat_wilson_vec( (su3_matrix *)(gen_pt[6][i]),
	    (wilson_vector *)(gen_pt[2][i]), &sittemp );
	mult_mat_wilson_vec( &(s->link[mu]), &sittemp, &ltemp );

	mult_mat_wilson_vec( (su3_matrix *)(gen_pt[1][i]),
	    (wilson_vector *)(gen_pt[3][i]), &rtemp );
	mult_adj_mat_wilson_vec( (su3_matrix *)(gen_pt[0][i]), &rtemp,
	    &sittemp);
	mult_sigma_mu_nu( &sittemp, &rtemp, mu, nu );

	su3_projector_w( &rtemp, &ltemp, &tmpmat );
	sub_su3_matrix( ((su3_matrix *)F_PT(s,mat)),
	    &tmpmat, ((su3_matrix *)F_PT(s,mat)) );
    }

    cleanup_gather(tag[0]);
    cleanup_gather(tag[1]);
    cleanup_general_gather(tag[2]);
    cleanup_general_gather(tag[3]);
    cleanup_gather(tag[4]);
    cleanup_gather(tag[5]);
    cleanup_general_gather(tag[6]);

} /* end udadu_mu_nu */
