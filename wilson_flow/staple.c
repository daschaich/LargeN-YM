// -----------------------------------------------------------------
// Construct staples
// Use tempmatf for temporary storage
#include "wflow_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Add staples defined by directions dir and dir2
void directional_staple(int dir, int dir2, field_offset lnk1,
                        field_offset lnk2, matrix_f *stp) {

  register int i;
  register site *s;
  msg_tag *tag0, *tag1, *tag2;
  matrix_f tmat;

  // Get blocked_link[dir2] from direction dir
  tag0 = start_gather_site(lnk2, sizeof(matrix_f), dir,
                           EVENANDODD, gen_pt[0]);

  // Get blocked_link[dir] from direction dir2
  tag1 = start_gather_site(lnk1, sizeof(matrix_f), dir2,
                           EVENANDODD, gen_pt[1]);

  // Start working on the lower staple while we wait for the gathers
  // The lower staple is prepared at x-dir2 and stored in tempmatf,
  // then gathered to x
  FORALLSITES(i, s) {
    mult_an_f((matrix_f *)F_PT(s, lnk2), (matrix_f *)F_PT(s, lnk1),
              &(tempmatf[i]));
  }

  // Finish lower staple
  wait_gather(tag0);
  wait_gather(tag1);
  FORALLSITES(i, s) {
    mult_nn_f(&(tempmatf[i]), (matrix_f *)gen_pt[0][i], &tmat);
    mat_copy_f(&tmat, &(tempmatf[i]));       // Overwrite tempmatf
  }

  // Gather staple from direction -dir2 to "home" site
  tag2 = start_gather_field(tempmatf, sizeof(matrix_f),
                            OPP_DIR(dir2), EVENANDODD, gen_pt[2]);

  // Calculate upper staple, add it
  FORALLSITES(i, s) {
    mult_nn_f((matrix_f *)F_PT(s, lnk2), (matrix_f *)gen_pt[1][i], &tmat);
    mult_na_sum_f(&tmat, (matrix_f *)gen_pt[0][i], &(stp[i]));
  }

  // Finally add the lower staple
  wait_gather(tag2);
  FORALLSITES(i, s)
    sum_mat_f((matrix_f *)gen_pt[2][i], &(stp[i]));

  cleanup_gather(tag0);
  cleanup_gather(tag1);
  cleanup_gather(tag2);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Sum staples for direction dir over all other directions
void staple(matrix_f *stp[NDIMS]) {
  register int i;
  register site *s;
  int dir, dir2;

  FORALLUPDIR(dir) {
    FORALLSITES(i, s)
      clear_mat_f(&(stp[dir][i]));

    FORALLUPDIR(dir2) {
      if (dir == dir2)
        continue;
      directional_staple(dir, dir2, F_OFFSET(linkf[dir]),
                         F_OFFSET(linkf[dir2]), stp[dir]);
    }
  }
}
// -----------------------------------------------------------------
