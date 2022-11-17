// -----------------------------------------------------------------
// Construct staples
// Use tempmat for temporary storage
#include "wflow_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Add staples defined by directions dir and dir2
void directional_staple(int dir, int dir2, field_offset lnk1,
                        field_offset lnk2, matrix *stp) {

  register int i;
  register site *s;
  msg_tag *tag, *tag2, *tag3;
  matrix tmat;

  // Get link[dir2] from direction dir
  tag = start_gather_site(lnk2, sizeof(matrix), dir,
                          EVENANDODD, gen_pt[0]);

  // Get link[dir] from direction dir2
  tag2 = start_gather_site(lnk1, sizeof(matrix), dir2,
                           EVENANDODD, gen_pt[1]);

  // Start working on the lower staple while we wait for the gathers
  // The lower staple is prepared at x-dir2 and stored in tempmat,
  // then gathered to x
  FORALLSITES(i, s)
    mult_an((matrix *)F_PT(s, lnk2), (matrix *)F_PT(s, lnk1), &(tempmat[i]));

  // Finish lower staple
  wait_gather(tag);
  wait_gather(tag2);
  FORALLSITES(i, s) {
    mult_nn(&(tempmat[i]), (matrix *)gen_pt[0][i], &tmat);
    mat_copy(&tmat, &(tempmat[i]));       // Overwrite tempmat
  }

  // Gather staple from direction -dir2 to "home" site
  tag3 = start_gather_field(tempmat, sizeof(matrix),
                            OPP_DIR(dir2), EVENANDODD, gen_pt[2]);

  // Calculate upper staple, add it
  FORALLSITES(i, s) {
    mult_nn((matrix *)F_PT(s, lnk2), (matrix *)gen_pt[1][i], &tmat);
    mult_na_sum(&tmat, (matrix *)gen_pt[0][i], &(stp[i]));
  }

  // Finally add the lower staple
  wait_gather(tag3);
  FORALLSITES(i, s)
    sum_mat((matrix *)gen_pt[2][i], &(stp[i]));

  cleanup_gather(tag);
  cleanup_gather(tag2);
  cleanup_gather(tag3);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Sum staples for direction dir over all other directions
void staple(matrix *stp[NDIMS]) {
  register int i;
  register site *s;
  int dir, dir2;

  FORALLUPDIR(dir) {
    FORALLSITES(i, s)
      clear_mat(&(stp[dir][i]));

    FORALLUPDIR(dir2) {
      if (dir == dir2)
        continue;
      directional_staple(dir, dir2, F_OFFSET(linkf[dir]),
                         F_OFFSET(linkf[dir2]), stp[dir]);
    }
  }
}
// -----------------------------------------------------------------
