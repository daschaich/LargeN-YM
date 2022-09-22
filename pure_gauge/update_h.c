// -----------------------------------------------------------------
// Update the momentum matrices
// Use tempmatf and tempmatf2 for temporary storage
#include "pg_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Update the momenta --- note the difference rather than sum
void update_anti_hermitian(site *s, int dir, Real eps, matrix_f *force) {
  matrix_f tmat;

  uncompress_anti_hermitian(&(s->mom[dir]), &tmat);
  scalar_mult_dif_mat_f(force, eps, &tmat);
  make_anti_hermitian(&tmat, &(s->mom[dir]));
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Update the momenta with the gauge force
double update_h(Real eps) {
  register int i, dir, dir2;
  register site *s;
  register Real ebN = eps * beta * one_ov_N;
  double norm = 0.0;
  msg_tag *tag0, *tag1, *tag2;
  int start;
  matrix_f tmat;

  // Loop over directions, update mom[dir]
  FORALLUPDIR(dir) {
    start = 1; // Indicates staple sum (in tempmat2) not initialized

    // Loop over other directions
    // Compute force from plaquettes in the dir, dir2 plane
    FORALLUPDIR(dir2) {
      if (dir2 == dir)
        continue;

      // Get linkf[dir2] from direction dir
      tag0 = start_gather_site(F_OFFSET(linkf[dir2]),
                               sizeof(matrix_f),
                               dir, EVENANDODD, gen_pt[0]);

      // Start gather for the "upper staple"
      tag2 = start_gather_site(F_OFFSET(linkf[dir]),
                               sizeof(matrix_f),
                               dir2, EVENANDODD, gen_pt[2]);

      // Begin the computation "at the dir2DOWN point"
      // We will later gather the intermediate result "to the home point"
      wait_gather(tag0);
      FORALLSITES(i, s) {
        mult_an_f(&(s->linkf[dir2]), &(s->linkf[dir]), &tmat);
        mult_nn_f(&tmat, (matrix_f *)gen_pt[0][i], &(tempmatf[i]));
      }

      // Gather lower staple "up to home site"
      tag1 = start_gather_field(tempmatf, sizeof(matrix_f),
                               OPP_DIR(dir2), EVENANDODD, gen_pt[1]);

      // The "upper" staple
      // One of the links has already been gathered,
      // since it was used in computing
      // the "lower" staple of the site above (in dir2)
      wait_gather(tag2);
      if (start) {  // Initialize staple sum in tempmatf2
        FORALLSITES(i, s) {
          mult_nn_f(&(s->linkf[dir2]), (matrix_f *)gen_pt[2][i], &tmat);
          mult_na_f(&tmat, (matrix_f *)gen_pt[0][i], &(tempmatf2[i]));
        }
        start = 0;
      }
      else {
        FORALLSITES(i, s) {
          mult_nn_f(&(s->linkf[dir2]), (matrix_f *)gen_pt[2][i], &tmat);
          mult_na_sum_f(&tmat, (matrix_f *)gen_pt[0][i], &(tempmatf2[i]));

        }
      }

      wait_gather(tag1);
      FORALLSITES(i, s)
        sum_mat_f((matrix_f *)gen_pt[1][i], &(tempmatf2[i]));
      cleanup_gather(tag0);
      cleanup_gather(tag1);
      cleanup_gather(tag2);
    }

    // Now multiply the staple sum by the link, then update momentum
    FORALLSITES(i, s) {
      mult_na_f(&(s->linkf[dir]), &(tempmatf2[i]), &tmat);
      update_anti_hermitian(s, dir, ebN, &tmat);
      realtrace_sum_f(&tmat, &tmat, &norm);
    }
  }
  g_doublesum(&norm);
  return (ebN * sqrt(norm) / (double)volume);
}

// Adds gaussian window term to force
// Should be possible to merge with routine above
double update_h_const(Real eps, double Eint, double a) {
#if 0
#ifdef LLR
  register int i, dir, dir2;
  register site *s;
  register Real ebN = eps * beta * one_ov_N;
  msg_tag *tag0, *tag1, *tag2;
  int start;
  matrix_f tmat, tmat2, tmat3, tmat4;
  double norm = 0.0;
  double CurrentEnergy = U_action();

  // Loop over directions, update mom[dir]
  FORALLUPDIR(dir) {
    start = 1; // Indicates staple sum (in tempmat2) not initialized

    // Loop over other directions
    // Compute force from plaquettes in the dir, dir2 plane
    FORALLUPDIR(dir2) {
      if (dir2 == dir)
        continue;

      // Get linkf[dir2] from direction dir
      tag0 = start_gather_site(F_OFFSET(linkf[dir2]),
                               sizeof(matrix_f),
                               dir, EVENANDODD, gen_pt[0]);

      // Start gather for the "upper staple"
      tag2 = start_gather_site(F_OFFSET(linkf[dir]),
                               sizeof(matrix_f),
                               dir2, EVENANDODD, gen_pt[2]);

      // Begin the computation "at the dir2DOWN point"
      // We will later gather the intermediate result "to the home point"
      wait_gather(tag0);
      FORALLSITES(i, s) {
        mult_an_f(&(s->linkf[dir2]), &(s->linkf[dir]), &tmat);
        mult_nn_f(&tmat, (matrix_f *)gen_pt[0][i], &(tempmatf[i]));
      }

      // Gather lower staple "up to home site"
      tag1 = start_gather_field(tempmatf, sizeof(matrix_f),
                               OPP_DIR(dir2), EVENANDODD, gen_pt[1]);

      // The "upper" staple
      // One of the links has already been gathered,
      // since it was used in computing
      // the "lower" staple of the site above (in dir2)
      wait_gather(tag2);
      if (start) {  // Initialize staple sum in tempmatf2
        FORALLSITES(i, s) {
          mult_nn_f(&(s->linkf[dir2]), (matrix_f *)gen_pt[2][i], &tmat);
          mult_na_f(&tmat, (matrix_f *)gen_pt[0][i], &(tempmatf2[i]));
        }
        start = 0;
      }
      else{
        FORALLSITES(i, s) {
          mult_nn_f(&(s->linkf[dir2]), (matrix_f *)gen_pt[2][i], &tmat);
          mult_na_f(&tmat, (matrix_f *)gen_pt[0][i], &tmat2);
          sum_mat_f(&tmat2, &(tempmatf2[i]));

        }
      }

      wait_gather(tag1);
      FORALLSITES(i, s)
        sum_mat_f((matrix_f *)gen_pt[1][i], &(tempmatf2[i]));
      cleanup_gather(tag0);
      cleanup_gather(tag1);
      cleanup_gather(tag2);
    }

    // Now multiply the staple sum by the link, then update momentum
    FORALLSITES(i, s) {
      mult_na_f(&(s->linkf[dir]), &(tempmatf2[i]), &tmat);
      uncompress_anti_hermitian(&(s->mom[dir]), &tmat2);
      scalar_mult_mat_f(&tmat,(CurrentEnergy-Eint-delta*0.5)/pow(delta,2),&tmat3);
      scalar_mult_add_mat_f(&tmat3, &tmat, a, &tmat4);
      scalar_mult_add_mat_f(&tmat2, &tmat4, -1.0 * ebN, &(tempmatf2[i]));
      make_anti_hermitian(&(tempmatf2[i]), &(s->mom[dir]));
      norm += (double)realtrace_f(&tmat, &tmat);
    }
  }

  g_doublesum(&norm);
  return ebN * sqrt(norm) / (double)volume;
#endif
#endif
  return(-99);
}
// -----------------------------------------------------------------
