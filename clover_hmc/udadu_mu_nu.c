/* This version uses gathers to get the neighbors */
/* Compute U dA/dU, part of the fermion force term, given mu & nu.

   The udadul terms are left-multiplied by the transpose of a
   wilson_vector, so we take the hermitian conjugate here, so we
   can operate on the vector from the left.  Essentially:
     x^* A = (A^\dagger x)^*
   The udadur terms are right-multiplied by a wilson_vector.

   The force term consists of only the right half of the "clover".
*/
#include "cl_dyn_includes.h"

// lsrc and rsrc are input
// Output is the contribution to iH_dot, put into tempmat
// Uses tempwvec; can't touch psi, chi, p and mp, which are busy
void udadu_mu_nu(wilson_vector *lsrc, wilson_vector *rsrc, int mu, int nu) {
  register int i;
  register site *s;
  int disp[4] = {0, 0, 0, 0};
  msg_tag *tag[8];
  wilson_vector ltemp, rtemp;   /* input to su3_projector_w() */
  wilson_vector sittemp;
  su3_matrix tmat;      /* for middle steps */

  /**********  First work on upper leaf  **********/
  /* get link[nu] from direction +mu */
  tag[0] = start_gather_site(F_OFFSET(link[nu]), sizeof(su3_matrix),
                             mu, EVENANDODD, gen_pt[0]);

  /* get link[mu] from direction +nu */
  tag[1] = start_gather_site(F_OFFSET(link[mu]), sizeof(su3_matrix),
                             nu, EVENANDODD, gen_pt[1]);

  /*****  Upper leaf - same parity *****/
  /*get sources from +mu +nu */
  disp[mu] = 1;
  disp[nu] = 1;
  tag[2] = start_general_gather_field(lsrc, sizeof(wilson_vector),
                                      disp, EVENANDODD, gen_pt[2]);

  wait_general_gather(tag[2]);
  tag[3] = start_general_gather_field(rsrc, sizeof(wilson_vector),
                                      disp, EVENANDODD, gen_pt[3]);

  /* ltemp = lsrc, rtemp = plaq*rsrc */
  wait_gather(tag[0]);
  wait_gather(tag[1]);
  FORALLSITES(i, s) {
    mult_adj_mat_wilson_vec(&(s->link[nu]), &(rsrc[i]), &rtemp);
    mult_adj_mat_wilson_vec((su3_matrix *)(gen_pt[1][i]), &rtemp, &sittemp);
    mult_mat_wilson_vec((su3_matrix *)(gen_pt[0][i]), &sittemp, &rtemp);
    mult_mat_wilson_vec(&(s->link[mu]), &rtemp, &sittemp);
    mult_sigma_mu_nu(&sittemp, &rtemp, mu, nu);

    // Initialize tempmat
    su3_projector_w(&rtemp, &(lsrc[i]), &(tempmat[i]));
  }

  /* ltemp = upperleftcorner*lsrc, rtemp = lowerrightcorner*rsrc */
  wait_general_gather(tag[3]);
  FORALLSITES(i, s) {
    mult_mat_wilson_vec((su3_matrix *)(gen_pt[1][i]),
                        (wilson_vector *)(gen_pt[2][i]), &sittemp);
    mult_mat_wilson_vec(&(s->link[nu]), &sittemp, &ltemp);

    mult_mat_wilson_vec((su3_matrix *)(gen_pt[0][i]),
                        (wilson_vector *)(gen_pt[3][i]), &rtemp);
    mult_mat_wilson_vec(&(s->link[mu]), &rtemp, &sittemp);
    mult_sigma_mu_nu(&sittemp, &rtemp, mu, nu);

    su3_projector_w(&rtemp, &ltemp, &tmat);
    sum_su3_matrix(&tmat, &(tempmat[i]));
  }
  cleanup_general_gather(tag[2]);
  cleanup_general_gather(tag[3]);

  /***** upper leaf - different parity *****/

  /* get sources from +nu */
  tag[4] = start_gather_field(lsrc, sizeof(wilson_vector),
                              nu, EVENANDODD, gen_pt[4]);
  tag[5] = start_gather_field(rsrc, sizeof(wilson_vector),
                              nu, EVENANDODD, gen_pt[5]);

  /* ltemp = link[nu]*lsrc, rtemp = staple*rsrc */
  wait_gather(tag[4]);
  wait_gather(tag[5]);
  FORALLSITES(i, s) {
    mult_mat_wilson_vec(&(s->link[nu]),
                        (wilson_vector *)(gen_pt[4][i]), &ltemp);

    mult_adj_mat_wilson_vec((su3_matrix *)(gen_pt[1][i]),
                            (wilson_vector *)(gen_pt[5][i]), &sittemp);
    mult_mat_wilson_vec((su3_matrix *)(gen_pt[0][i]), &sittemp, &rtemp);
    mult_mat_wilson_vec(&(s->link[mu]), &rtemp, &sittemp);
    mult_sigma_mu_nu(&sittemp, &rtemp, mu, nu);
    su3_projector_w(&rtemp, &ltemp, &tmat);
    sum_su3_matrix(&tmat, &(tempmat[i]));
  }
  cleanup_gather(tag[4]);
  cleanup_gather(tag[5]);

  /* get sources from +mu */
  tag[4] = start_gather_field(lsrc, sizeof(wilson_vector),
                              mu, EVENANDODD, gen_pt[4]);
  tag[5] = start_gather_field(rsrc, sizeof(wilson_vector),
                              mu, EVENANDODD, gen_pt[5]);

  /* ltemp = staple*lsrc, rtemp = link[mu]*rsrc */
  wait_gather(tag[4]);
  wait_gather(tag[5]);
  FORALLSITES(i, s) {
    mult_adj_mat_wilson_vec((su3_matrix *)(gen_pt[0][i]),
                            (wilson_vector *)(gen_pt[4][i]), &ltemp);
    mult_mat_wilson_vec((su3_matrix *)(gen_pt[1][i]), &ltemp, &sittemp);
    mult_mat_wilson_vec(&(s->link[nu]), &sittemp, &ltemp);

    mult_mat_wilson_vec(&(s->link[mu]), (wilson_vector *)(gen_pt[5][i]),
                        &sittemp);
    mult_sigma_mu_nu(&sittemp, &rtemp, mu, nu);

    su3_projector_w(&rtemp, &ltemp, &tmat);
    sum_su3_matrix(&tmat, &(tempmat[i]));
  }
  cleanup_gather(tag[0]);
  cleanup_gather(tag[1]);

  /* Hold on to these sources since they appear in lower leaf below */

  /**********  The Lower leaf  **********/

  /* get link[nu] from direction -nu */
  tag[0] = start_gather_site(F_OFFSET(link[nu]), sizeof(su3_matrix),
                             OPP_DIR(nu), EVENANDODD, gen_pt[0]);

  /* get link[mu] from direction -nu */
  tag[1] = start_gather_site(F_OFFSET(link[mu]), sizeof(su3_matrix),
                             OPP_DIR(nu), EVENANDODD, gen_pt[1]);

  /* get link[nu] from direction -nu +mu; disp[mu] = 1 already */
  disp[nu] = -1;
  tag[6] = start_general_gather_site(F_OFFSET(link[nu]), sizeof(su3_matrix),
                                     disp, EVENANDODD, gen_pt[6]);

  wait_gather(tag[0]);
  wait_gather(tag[1]);
  wait_general_gather(tag[6]);

  /***** lower leaf, different parity *****/

  /* ltemp = link[mu]*lsrc, rtemp = staple*rsrc */
  FORALLSITES(i, s) {
    mult_mat_wilson_vec(&(s->link[mu]), (wilson_vector *)(gen_pt[4][i]),
                        &ltemp);

    mult_mat_wilson_vec((su3_matrix *)(gen_pt[6][i]),
                        (wilson_vector *)(gen_pt[5][i]), &sittemp);
    mult_mat_wilson_vec((su3_matrix *)(gen_pt[1][i]), &sittemp, &rtemp);
    mult_adj_mat_wilson_vec((su3_matrix *)(gen_pt[0][i]), &rtemp, &sittemp);
    mult_sigma_mu_nu(&sittemp, &rtemp, mu, nu);

    su3_projector_w(&rtemp, &ltemp, &tmat);
    dif_su3_matrix(&tmat, &(tempmat[i]));
  }
  cleanup_gather(tag[4]);
  cleanup_gather(tag[5]);

  /* get sources from -nu */
  tag[4] = start_gather_field(lsrc, sizeof(wilson_vector),
                              OPP_DIR(nu), EVENANDODD, gen_pt[4]);
  tag[5] = start_gather_field(rsrc, sizeof(wilson_vector),
                              OPP_DIR(nu), EVENANDODD, gen_pt[5]);

  /* ltemp = staple*lsrc, rtemp = link[nu]^dagger*rsrc */
  wait_gather(tag[4]);
  wait_gather(tag[5]);
  FORALLSITES(i, s) {
    mult_adj_mat_wilson_vec((su3_matrix *)(gen_pt[1][i]),
                            (wilson_vector *)(gen_pt[4][i]), &ltemp);
    mult_adj_mat_wilson_vec((su3_matrix *)(gen_pt[6][i]), &ltemp, &sittemp);
    mult_mat_wilson_vec(&(s->link[mu]), &sittemp, &ltemp);

    mult_adj_mat_wilson_vec((su3_matrix *)(gen_pt[0][i]),
                            (wilson_vector *)(gen_pt[5][i]), &sittemp);
    mult_sigma_mu_nu(&sittemp, &rtemp, mu, nu);

    su3_projector_w(&rtemp, &ltemp, &tmat);
    dif_su3_matrix(&tmat, &(tempmat[i]));
  }
  cleanup_gather(tag[4]);
  cleanup_gather(tag[5]);

  /***** lower leaf, same parity *****/

  /* get sources from +mu -nu */
  /* disp[mu] = 1;  disp[nu] = -1; already */
  tag[2] = start_general_gather_field(lsrc, sizeof(wilson_vector),
                                      disp, EVENANDODD, gen_pt[2]);
  wait_general_gather(tag[2]);
  tag[3] = start_general_gather_field(rsrc, sizeof(wilson_vector),
                                      disp, EVENANDODD, gen_pt[3]);

  /* ltemp = plaq*lsrc, rtemp = rsrc */
  FORALLSITES(i, s) {
    mult_mat_wilson_vec((su3_matrix *)(gen_pt[0][i]), &(lsrc[i]), &sittemp);
    mult_adj_mat_wilson_vec((su3_matrix *)(gen_pt[1][i]), &sittemp, &ltemp);
    mult_adj_mat_wilson_vec((su3_matrix *)(gen_pt[6][i]), &ltemp, &sittemp);
    mult_mat_wilson_vec(&(s->link[mu]), &sittemp, &ltemp);

    mult_sigma_mu_nu(&(rsrc[i]), &rtemp, mu, nu);

    su3_projector_w(&rtemp, &ltemp, &tmat);
    dif_su3_matrix(&tmat, &(tempmat[i]));
  }

  /* ltemp = upperrightcorner*lsrc, rtemp = lowerleftcorner*rsrc */
  wait_general_gather(tag[3]);
  FORALLSITES(i, s) {
    mult_adj_mat_wilson_vec((su3_matrix *)(gen_pt[6][i]),
                            (wilson_vector *)(gen_pt[2][i]), &sittemp);
    mult_mat_wilson_vec(&(s->link[mu]), &sittemp, &ltemp);

    mult_mat_wilson_vec((su3_matrix *)(gen_pt[1][i]),
                        (wilson_vector *)(gen_pt[3][i]), &rtemp);
    mult_adj_mat_wilson_vec((su3_matrix *)(gen_pt[0][i]), &rtemp, &sittemp);
    mult_sigma_mu_nu(&sittemp, &rtemp, mu, nu);

    su3_projector_w(&rtemp, &ltemp, &tmat);
    dif_su3_matrix(&tmat, &(tempmat[i]));
  }
  cleanup_gather(tag[0]);
  cleanup_gather(tag[1]);
  cleanup_general_gather(tag[2]);
  cleanup_general_gather(tag[3]);
  cleanup_general_gather(tag[6]);
}
// -----------------------------------------------------------------
