// -----------------------------------------------------------------
// Update the momentum matrices
// Use tempmat and tempmat2 for temporary storage
#include "pg_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Update the momenta --- note the difference rather than sum
void update_anti_hermitian(site *s, int dir, Real eps, matrix *force) {
  matrix tmat;

  uncompress_anti_hermitian(&(s->mom[dir]), &tmat);
  scalar_mult_dif_mat(force, eps, &tmat);
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
  matrix tmat;

  // Loop over directions, update mom[dir]
  FORALLUPDIR(dir) {
    start = 1; // Indicates staple sum (in tempmat2) not initialized

    // Loop over other directions
    // Compute force from plaquettes in the dir, dir2 plane
    FORALLUPDIR(dir2) {
      if (dir2 == dir)
        continue;

      // Get link[dir2] from direction dir
      tag0 = start_gather_site(F_OFFSET(link[dir2]), sizeof(matrix),
                               dir, EVENANDODD, gen_pt[0]);

      // Start gather for the "upper staple"
      tag2 = start_gather_site(F_OFFSET(link[dir]), sizeof(matrix),
                               dir2, EVENANDODD, gen_pt[2]);

      // Begin the computation "at the dir2DOWN point"
      // We will later gather the intermediate result "to the home point"
      wait_gather(tag0);
      FORALLSITES(i, s) {
        mult_an(&(s->link[dir2]), &(s->link[dir]), &tmat);
        mult_nn(&tmat, (matrix *)gen_pt[0][i], &(tempmat[i]));
      }

      // Gather lower staple "up to home site"
      tag1 = start_gather_field(tempmat, sizeof(matrix),
                               OPP_DIR(dir2), EVENANDODD, gen_pt[1]);

      // The "upper" staple
      // One of the links has already been gathered,
      // since it was used in computing
      // the "lower" staple of the site above (in dir2)
      wait_gather(tag2);
      if (start) {  // Initialize staple sum in tempmat2
        FORALLSITES(i, s) {
          mult_nn(&(s->link[dir2]), (matrix *)gen_pt[2][i], &tmat);
          mult_na(&tmat, (matrix *)gen_pt[0][i], &(tempmat2[i]));
        }
        start = 0;
      }
      else {
        FORALLSITES(i, s) {
          mult_nn(&(s->link[dir2]), (matrix *)gen_pt[2][i], &tmat);
          mult_na_sum(&tmat, (matrix *)gen_pt[0][i], &(tempmat2[i]));

        }
      }

      wait_gather(tag1);
      FORALLSITES(i, s)
        sum_mat((matrix *)gen_pt[1][i], &(tempmat2[i]));
      cleanup_gather(tag0);
      cleanup_gather(tag1);
      cleanup_gather(tag2);
    }

    // Now multiply the staple sum by the link, then update momentum
    FORALLSITES(i, s) {
      mult_na(&(s->link[dir]), &(tempmat2[i]), &tmat);
      update_anti_hermitian(s, dir, ebN, &tmat);
      realtrace_sum(&tmat, &tmat, &norm);
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
  matrix tmat, tmat2, tmat3, tmat4;
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

      // Get link[dir2] from direction dir
      tag0 = start_gather_site(F_OFFSET(link[dir2]), sizeof(matrix),
                               dir, EVENANDODD, gen_pt[0]);

      // Start gather for the "upper staple"
      tag2 = start_gather_site(F_OFFSET(link[dir]), sizeof(matrix),
                               dir2, EVENANDODD, gen_pt[2]);

      // Begin the computation "at the dir2DOWN point"
      // We will later gather the intermediate result "to the home point"
      wait_gather(tag0);
      FORALLSITES(i, s) {
        mult_an(&(s->link[dir2]), &(s->link[dir]), &tmat);
        mult_nn(&tmat, (matrix *)gen_pt[0][i], &(tempmat[i]));
      }

      // Gather lower staple "up to home site"
      tag1 = start_gather_field(tempmat, sizeof(matrix),
                               OPP_DIR(dir2), EVENANDODD, gen_pt[1]);

      // The "upper" staple
      // One of the links has already been gathered,
      // since it was used in computing
      // the "lower" staple of the site above (in dir2)
      wait_gather(tag2);
      if (start) {  // Initialize staple sum in tempmat2
        FORALLSITES(i, s) {
          mult_nn(&(s->link[dir2]), (matrix *)gen_pt[2][i], &tmat);
          mult_na(&tmat, (matrix *)gen_pt[0][i], &(tempmat2[i]));
        }
        start = 0;
      }
      else{
        FORALLSITES(i, s) {
          mult_nn(&(s->link[dir2]), (matrix *)gen_pt[2][i], &tmat);
          mult_na(&tmat, (matrix *)gen_pt[0][i], &tmat2);
          sum_mat(&tmat2, &(tempmat2[i]));

        }
      }

      wait_gather(tag1);
      FORALLSITES(i, s)
        sum_mat((matrix *)gen_pt[1][i], &(tempmat2[i]));
      cleanup_gather(tag0);
      cleanup_gather(tag1);
      cleanup_gather(tag2);
    }

    // Now multiply the staple sum by the link, then update momentum
    FORALLSITES(i, s) {
      mult_na(&(s->link[dir]), &(tempmat2[i]), &tmat);
      uncompress_anti_hermitian(&(s->mom[dir]), &tmat2);
      scalar_mult_mat(&tmat,(CurrentEnergy-Eint-delta*0.5)/pow(delta,2),&tmat3);
      scalar_mult_add_mat(&tmat3, &tmat, a, &tmat4);
      scalar_mult_add_mat(&tmat2, &tmat4, -1.0 * ebN, &(tempmat2[i]));
      make_anti_hermitian(&(tempmat2[i]), &(s->mom[dir]));
      norm += (double)realtrace(&tmat, &tmat);
    }
  }

  g_doublesum(&norm);
  return ebN * sqrt(norm) / (double)volume;
#endif
#endif
  return(-99);
}
// -----------------------------------------------------------------
