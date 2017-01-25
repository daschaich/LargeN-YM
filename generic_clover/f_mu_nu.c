// -----------------------------------------------------------------
/* This version uses (general) gathers to get the neighbors */

// Compute F_{\mu\nu} used in the clover fermion action

/*  f_mn[site] = (G_mn - G_mn^dag) / 8  (antihermitian) where

    G_mn is given by the product of link matrices in the pattern shown.

          <-----^                  <-----^
          |     |                  |     |
          |     |                  |     |
  G_mn =  O----->   +  <-----O  +  v---->O  +  O<----^
                       |     |                 |     |
                       |     |                 |     |
     nu  ^             v----->                 v----->
         |
         O---> mu

    The "site" is indicated by O.

    If we write U_mu = exp[i A_mu] then

    f_mn \approx i F_{\mu\nu}
 */

#include "generic_clover_includes.h"

// Use tempmat, tempmat2 and staple for temporary storage
void f_mu_nu(matrix f_mn[], int mu, int nu) {
  register int i;
  register site *s;
  int disp[4] = {0, 0, 0, 0};  /* displacement for general gather */
  int order_flag;
  matrix tmat;
  msg_tag *tag0, *tag1, *tag2, *tag3, *tag4;

  /* Want mu < nu, so that only nu can be TUP! */
  if (mu > nu) {
    i = mu;
    mu = nu;
    nu = i;
    order_flag = 1;
  }
  else
    order_flag = 0;

  /* get link[nu] from direction +mu */
  tag0 = start_gather_site(F_OFFSET(link[nu]), sizeof(matrix),
                           mu, EVENANDODD, gen_pt[0]);

  /* get link[mu] from direction +nu */
  tag1 = start_gather_site(F_OFFSET(link[mu]), sizeof(matrix),
                           nu, EVENANDODD, gen_pt[1]);

  /* Make one corner with link[nu]^dag link[mu] */
  FORALLSITES(i, s)
    mult_an(&(s->link[nu]), &(s->link[mu]), &(tempmat[i]));

  /* Make one corner with link[nu](x+mu) link[mu](x+nu)^dagger
     and multiply the two corners together in the two different ways */
  /* Note f_mn is here used as a temporary! */
  wait_gather(tag0);
  wait_gather(tag1);
  FORALLSITES(i, s) {
    mult_na((matrix *)(gen_pt[0][i]), (matrix *)gen_pt[1][i],
                &(f_mn[i]));
    mult_nn(&(tempmat[i]), &(f_mn[i]), &(tempmat2[i]));
    mult_nn(&(f_mn[i]), &(tempmat[i]), &(staple[i]));
  }

  /* tempmat2 is the plaquette +mu -nu and must be gathered from -nu */
  tag2 = start_gather_field(tempmat2, sizeof(matrix),
                            OPP_DIR(nu), EVENANDODD, gen_pt[2]);

  /* staple is the plaquette -mu +nu and must be gather from -mu */
  tag3 = start_gather_field(staple, sizeof(matrix),
                            OPP_DIR(mu), EVENANDODD, gen_pt[3]);

  /* Now make +mu +nu plaquette and put in f_mn */
  FORALLSITES(i, s) {
    mult_nn(&(s->link[mu]), &(f_mn[i]), &tmat);
    mult_na(&tmat, &(s->link[nu]), &(f_mn[i]));
  }

  /* Now gather +mu -nu plaquette and add to f_mn */
  wait_gather(tag2);
  FORALLSITES(i, s)
    sum_mat((matrix *)(gen_pt[2][i]), &(f_mn[i]));

  cleanup_gather(tag2);
  FORALLSITES(i, s) {
    mult_an((matrix *)(gen_pt[1][i]), &(tempmat[i]), &tmat);
    mult_nn(&tmat, (matrix *)(gen_pt[0][i]), &(tempmat2[i]));
  }
  cleanup_gather(tag0);
  cleanup_gather(tag1);

  /* tempmat2 is now plaquette -mu -nu and must be gathered with
     displacement -mu-nu */
  disp[mu] = -1;
  disp[nu] = -1;
  tag4 = start_general_gather_field(tempmat2, sizeof(matrix),
                                    disp, EVENANDODD, gen_pt[4]);

  /* Now gather -mu +nu plaquette and add to f_mn */
  wait_gather(tag3);
  FORALLSITES(i, s)
    sum_mat((matrix *)(gen_pt[3][i]), &(f_mn[i]));
  cleanup_gather(tag3);

  /* Finally gather -mu -nu plaquette and add to f_mn */
  wait_general_gather(tag4);
  FORALLSITES(i, s)
    sum_mat((matrix *)(gen_pt[4][i]), &(f_mn[i]));
  cleanup_general_gather(tag4);

  // Final factor of 1 / 8 on (f_mn - f_mn^dag)
  FORALLSITES(i, s) {
    adjoint(&(f_mn[i]), &tmat);
    if (order_flag == 0)
      sub_mat(&(f_mn[i]), &tmat, &tmat);
    else
      dif_mat(&(f_mn[i]), &tmat);
    scalar_mult_mat(&tmat, 0.125, &(f_mn[i]));
  }
}
// -----------------------------------------------------------------
