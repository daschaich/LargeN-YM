// -----------------------------------------------------------------
// Compute the staple
#include "pg_includes.h"
#include "../include/loopend.h"

void dsdu_qhb(int dir, int parity) {
  register int i, dir2, otherparity = 0;
  register site *s;
  msg_tag *tag0, *tag1, *tag2, *tag3;
  int start;
  matrix tmat;

  switch(parity) {
    case EVEN:       otherparity = ODD;        break;
    case ODD:        otherparity = EVEN;       break;
    case EVENANDODD: otherparity = EVENANDODD; break;
  }

  // Loop over other directions,
  // computing force from plaquettes in the dir,dir2 plane
  start = 1; // Indicates staple sum not initialized
  FORALLUPDIR(dir2) {
    if (dir2 == dir)
      continue;

    // Get link[dir2] from direction dir on other parity
    tag0 = start_gather_site(F_OFFSET(link[dir2]), sizeof(matrix),
                             dir, otherparity, gen_pt[0]);

    // Get link[dir2] from direction dir
    tag1 = start_gather_site(F_OFFSET(link[dir2]), sizeof(matrix),
                             dir, parity, gen_pt[1]);

    // Get link[dir] from direction dir2
    tag2 = start_gather_site(F_OFFSET(link[dir]), sizeof(matrix),
                             dir2, parity, gen_pt[2]);

    // Lower staple (computed at backward site)
    wait_gather(tag0);
    FORSOMEPARITY(i, s, otherparity) {
      mult_an(&(s->link[dir2]), &(s->link[dir]), &tmat);
      mult_nn(&tmat, (matrix *)gen_pt[0][i], &(s->tempmat));
    } END_LOOP
    cleanup_gather(tag0);
    // Get tempmat from direction -dir2
    tag3 = start_gather_site(F_OFFSET(tempmat), sizeof(matrix),
                             OPP_DIR(dir2), parity, gen_pt[3]);

    // Upper staple
    wait_gather(tag1);
    wait_gather(tag2);
    if (start) {  // This is the first contribution to staple
      FORSOMEPARITY(i, s, parity) {
        mult_nn(&(s->link[dir2]), (matrix *)gen_pt[2][i], &tmat);
        mult_na(&tmat, (matrix *)gen_pt[1][i], &(s->staple));
      } END_LOOP
      start = 0;
    }
    else {
      FORSOMEPARITY(i, s, parity) {
        mult_nn(&(s->link[dir2]), (matrix *)gen_pt[2][i], &tmat);
        mult_na_sum(&tmat, (matrix *)gen_pt[1][i], &(s->staple));
      } END_LOOP
    }
    cleanup_gather(tag1);
    cleanup_gather(tag2);

    // Add lower staple
    wait_gather(tag3);
    FORSOMEPARITY(i, s, parity) {
      sum_mat((matrix *)gen_pt[3][i], &(s->staple));
    } END_LOOP
    cleanup_gather(tag3);
  }
}
// -----------------------------------------------------------------
