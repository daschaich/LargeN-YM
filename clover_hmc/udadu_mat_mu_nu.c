/* calculates udadu_mu_nu and multiplies by an matrix
   for the Tr_{dirac} \sigma_mu_nu A_odd^{-1} U dA_odd/dU term of
   the fermion force
*/
// Must be called after tr_sigma_ldu_mu_nu constructs tempmat
// Put answer into tempmat2

#include "cl_dyn_includes.h"

void udadu_mat_mu_nu(int mu, int nu) {
  int disp[4];
  register int i;
  register site *s;
  msg_tag *tag[8];
  matrix tmat, utmp, uAcntr;    /* for middle steps */

  /**********  First work on upper leaf  **********/
  /* get link[nu] from direction +mu */
  tag[0] = start_gather_site(F_OFFSET(link[nu]), sizeof(matrix),
                             mu, EVENANDODD, gen_pt[0]);

  /* get link[mu] from direction +nu */
  tag[1] = start_gather_site(F_OFFSET(link[mu]), sizeof(matrix),
                             nu, EVENANDODD, gen_pt[1]);

  /*****  Upper leaf - same parity *****/
  /* get tempmat from +mu +nu */
  FORALLUPDIR(i)
    disp[i]=0;
  disp[mu] = 1;
  disp[nu] = 1;
  tag[2] = start_general_gather_field(tempmat, sizeof(matrix),
                                      disp, ODD, gen_pt[2]);

  wait_gather(tag[0]);
  wait_gather(tag[1]);

  /* utmp = plaquette */
  // Initialize tempmat2 on odd sites
  FORODDSITES(i, s) {
    mult_nn(&(s->link[mu]), (matrix *)(gen_pt[0][i]), &utmp);
    mult_na(&utmp, (matrix *)(gen_pt[1][i]), &tmat);
    mult_na(&tmat, &(s->link[nu]), &utmp);
    mult_nn(&utmp, &(tempmat[i]), &(tempmat2[i]));
  }

  wait_general_gather(tag[2]);

  /* utmp = upperleftcorner, uAcntr = lowerrightcorner*tempmat */
  FORODDSITES(i, s) {
    mult_nn(&(s->link[mu]), (matrix *)(gen_pt[0][i]), &tmat);
    mult_nn(&tmat, (matrix *)(gen_pt[2][i]), &uAcntr);
    mult_nn(&(s->link[nu]), (matrix *)(gen_pt[1][i]), &utmp);
    mult_na_sum(&uAcntr, &utmp, &(tempmat2[i]));
  }

  /***** upper leaf - different parity *****/

  /* get tempmat from +nu */
  tag[3] = start_gather_field(tempmat, sizeof(matrix), nu, EVEN, gen_pt[3]);
  wait_gather(tag[3]);

  /* uAcntr = staple*tempmat, utmp = link[nu]^dagger */
  // Initialize tempmat2 on even sites
  FOREVENSITES(i, s) {
    mult_nn(&(s->link[mu]), (matrix *)(gen_pt[0][i]), &tmat);
    mult_na(&tmat, (matrix *)(gen_pt[1][i]), &utmp);
    mult_nn(&utmp, (matrix *)(gen_pt[3][i]), &uAcntr);
    mult_na(&uAcntr, &(s->link[nu]), &(tempmat2[i]));
  }

  cleanup_gather(tag[3]);
  /* get tempmat from +mu */
  tag[3] = start_gather_field(tempmat, sizeof(matrix), mu, EVEN, gen_pt[3]);
  wait_gather(tag[3]);

  /* uAcntr = link[mu]*tempmat, utmp = staple */
  FOREVENSITES(i, s) {
    mult_nn(&(s->link[mu]), (matrix *)(gen_pt[3][i]), &uAcntr);
    mult_nn(&(s->link[nu]), (matrix *)(gen_pt[1][i]), &tmat);
    mult_na(&tmat, (matrix *)(gen_pt[0][i]), &utmp);
    mult_na_sum(&uAcntr, &utmp, &(tempmat2[i]));
  }

  /**********  The Lower leaf  **********/

  cleanup_gather(tag[0]);
  cleanup_gather(tag[1]);
  /* get link[nu] from direction -nu */
  tag[0] = start_gather_site(F_OFFSET(link[nu]), sizeof(matrix),
                             OPP_DIR(nu), EVENANDODD, gen_pt[0]);

  /* get link[mu] from direction -nu */
  tag[1] = start_gather_site(F_OFFSET(link[mu]), sizeof(matrix),
                                      OPP_DIR(nu), EVENANDODD, gen_pt[1]);

  /* get link[nu] from direction -nu +mu */
  /* disp[mu] = 1; already */
  disp[nu] = -1;
  tag[4] = start_general_gather_site(F_OFFSET(link[nu]), sizeof(matrix),
                                     disp, EVENANDODD, gen_pt[4]);

  wait_gather(tag[0]);
  wait_gather(tag[1]);
  wait_general_gather(tag[4]);

  /***** lower leaf, different parity *****/

  /* uAcntr = staple*tempmat, utmp = link[mu]^dagger */
  FOREVENSITES(i, s) {
    mult_an((matrix *)(gen_pt[0][i]), (matrix *)(gen_pt[1][i]),
                &uAcntr);
    mult_nn(&uAcntr, (matrix *)(gen_pt[4][i]), &tmat);
    mult_nn(&tmat, (matrix *)(gen_pt[3][i]), &uAcntr);
    mult_na_dif(&uAcntr, &(s->link[mu]), &(tempmat2[i]));
  }

  cleanup_gather(tag[3]);
  /* get tempmat from -nu */
  tag[3] = start_gather_field(tempmat, sizeof(matrix), OPP_DIR(nu),
                              EVEN, gen_pt[3]);
  wait_gather(tag[3]);

  /* uAcntr = link[nu]^dagger*tempmat, utmp=staple */
  FOREVENSITES(i, s) {
    mult_an((matrix *)(gen_pt[0][i]), (matrix *)(gen_pt[3][i]),
                &uAcntr);
    mult_nn((matrix *)(gen_pt[1][i]), (matrix *)(gen_pt[4][i]),
                &tmat);
    mult_na(&tmat, &(s->link[mu]), &utmp);
    mult_nn_dif(&uAcntr, &utmp, &(tempmat2[i]));
  }

  /***** lower leaf, same parity *****/

  cleanup_general_gather(tag[2]);
  /* get tempmat from +mu -nu */
  /* disp[mu] = 1;  disp[nu] = -1; already */
  tag[2] = start_general_gather_field(tempmat, sizeof(matrix), disp,
                                      ODD, gen_pt[2]);

  /* utmp= plaq */
  FORODDSITES(i, s) {
    mult_an((matrix *)(gen_pt[0][i]), (matrix *)(gen_pt[1][i]),
                &utmp);
    mult_nn(&utmp, (matrix *)(gen_pt[4][i]), &tmat);
    mult_na(&tmat, &(s->link[mu]), &utmp);

    mult_nn_dif(&(tempmat[i]), &utmp, &(tempmat2[i]));
  }

  wait_general_gather(tag[2]);

  /* uAcntr = lowerleftcorner*tempmat, utmp=upperrightcorner */
  FORODDSITES(i, s) {
    mult_an((matrix *)(gen_pt[0][i]), (matrix *)(gen_pt[1][i]),
                &tmat);
    mult_nn(&tmat, (matrix *)(gen_pt[2][i]), &uAcntr);
    mult_na(&(s->link[mu]), (matrix *)(gen_pt[4][i]), &utmp);
    mult_na_dif(&uAcntr, &utmp, &(tempmat2[i]));
  }

  cleanup_gather(tag[0]);
  cleanup_gather(tag[1]);
  cleanup_general_gather(tag[2]);
  cleanup_gather(tag[3]);
  cleanup_general_gather(tag[4]);
}
// -----------------------------------------------------------------
