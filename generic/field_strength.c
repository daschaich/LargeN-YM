// -----------------------------------------------------------------
// Compute SU(NCOL) field strength tensor
// Need to define six field strength indices
// Use tempmatf and tempmatf2 for temporary storage

// This cartoon shows how the plaquettes are calculated
// The path begins and ends at the 'O' at the corner
// F_munu is the sum of these plaquettes minus their adjoints

//  ^         --------<--------       --------<--------
//  |dir2     |               |       |               |
//  |         |               |       |               |
//  |         |               |       |               |
//  |         |               ^       |               ^
//  ------>   |               |       |               |
//    dir     |               |       |               |
//            |               |       |               |
//            |               |       |               |
//            -------->-------O       O------->--------
//
//            --------<-------O       O-------<--------
//            |               |       |               |
//            |               |       |               |
//            |               |       |               |
//            |               ^       |               ^
//            |               |       |               |
//            |               |       |               |
//            |               |       |               |
//            |               |       |               |
//            -------->--------       -------->--------

// Convention: try to use gen_pt[0] and mtag for links in direction dir
// These are gathered from +/- dir2

#include "generic_includes.h"

#define LINK_OFFSET(dir) link_src + sizeof(matrix_f) * (dir)
#define LINK(dir) (((matrix_f *)F_PT(s, link_src))[dir])
#define FS(component) (((matrix_f *)F_PT(s, field_dest))[component])

// link_src is offset for matrix_f linkf[4] in site struct
// field_dest is offset for matrix_f fieldstrength[6] in site struct
void make_field_strength(field_offset link_src, field_offset field_dest) {
  register int i, component, dir = -99, dir2 = -99;
  register site *s;
  int j;
  complex tc;
  matrix_f tmat, tmat2;
  msg_tag *mtag, *mtag2;

  for (component = FS_XY; component <= FS_ZT; component++) {
    switch(component) {
      case FS_XY: dir = XUP; dir2 = YUP; break;
      case FS_XZ: dir = XUP; dir2 = ZUP; break;
      case FS_YZ: dir = YUP; dir2 = ZUP; break;
      case FS_XT: dir = XUP; dir2 = TUP; break;
      case FS_YT: dir = YUP; dir2 = TUP; break;
      case FS_ZT: dir = ZUP; dir2 = TUP; break;
    }

    // +dir +dir2 plaquette
    mtag  = start_gather_site(LINK_OFFSET(dir), sizeof(matrix_f),
                              dir2, EVENANDODD, gen_pt[0]);
    mtag2 = start_gather_site(LINK_OFFSET(dir2), sizeof(matrix_f),
                              dir, EVENANDODD, gen_pt[1]);

    wait_gather(mtag);
    wait_gather(mtag2);
    FORALLSITES(i, s) {
      mult_nn_f(&LINK(dir), (matrix_f *)(gen_pt[1][i]), &tmat);
      mult_na_f(&tmat, (matrix_f *)(gen_pt[0][i]), &tmat2);
      mult_na_f(&tmat2, &LINK(dir2), &tmat);
      adjoint_f(&tmat, &tmat2);
      sub_mat_f(&tmat, &tmat2, &FS(component));
    }
    cleanup_gather(mtag2);

    // -dir +dir2 plaquette
    // Reuse link[dir] gather from dir2 corresponding to mtag
    FORALLSITES(i, s) {
      mult_an_f(&LINK(dir2), &LINK(dir), &tmat);
      mult_an_f((matrix_f *)(gen_pt[0][i]), &tmat, &(tempmatf[i]));
    }
    mtag2 = start_gather_field(tempmatf, sizeof(matrix_f),
                               OPP_DIR(dir), EVENANDODD, gen_pt[1]);

    wait_gather(mtag2);
    FORALLSITES(i, s) {
      mult_nn_f(&LINK(dir2), (matrix_f *)(gen_pt[1][i]), &tmat);
      adjoint_f(&tmat, &tmat2);
      sum_mat_f(&tmat, &FS(component));
      dif_mat_f(&tmat2, &FS(component));
    }
    cleanup_gather(mtag);
    cleanup_gather(mtag2);

    // -dir -dir2 plaquette
    mtag = start_gather_site(LINK_OFFSET(dir), sizeof(matrix_f),
                             OPP_DIR(dir), EVENANDODD, gen_pt[0]);
    mtag2 = start_gather_site(LINK_OFFSET(dir2), sizeof(matrix_f),
                              OPP_DIR(dir2), EVENANDODD, gen_pt[1]);

    wait_gather(mtag);
    wait_gather(mtag2);
    FORALLSITES(i, s) {
      mult_nn_f((matrix_f *)(gen_pt[0][i]), &LINK(dir2), &(tempmatf[i]));
      mult_nn_f((matrix_f *)(gen_pt[1][i]), &LINK(dir), &(tempmatf2[i]));
    }
    cleanup_gather(mtag);
    cleanup_gather(mtag2);

    mtag = start_gather_field(tempmatf, sizeof(matrix_f),
                              OPP_DIR(dir2), EVENANDODD, gen_pt[0]);
    mtag2 = start_gather_field(tempmatf2, sizeof(matrix_f),
                               OPP_DIR(dir), EVENANDODD, gen_pt[1]);

    wait_gather(mtag);
    wait_gather(mtag2);
    FORALLSITES(i, s) {
      mult_an_f((matrix_f *)(gen_pt[1][i]), (matrix_f *)(gen_pt[0][i]), &tmat);
      adjoint_f(&tmat, &tmat2);
      sum_mat_f(&tmat, &FS(component));
      dif_mat_f(&tmat2, &FS(component));
    }
    cleanup_gather(mtag);
    cleanup_gather(mtag2);

    // +dir -dir2 plaquette
    mtag2 = start_gather_site(LINK_OFFSET(dir2), sizeof(matrix_f),
                              dir, EVENANDODD, gen_pt[1]);

    wait_gather(mtag2);
    FORALLSITES(i, s) {
      mult_an_f(&LINK(dir2), &LINK(dir), &tmat);
      mult_nn_f(&tmat, (matrix_f *)(gen_pt[1][i]), &tempmatf[i]);
    }
    cleanup_gather(mtag2);

    mtag = start_gather_field(tempmatf, sizeof(matrix_f),
                              OPP_DIR(dir2), EVENANDODD, gen_pt[0]);
    wait_gather(mtag);
    FORALLSITES(i, s) {
      mult_na_f((matrix_f *)(gen_pt[0][i]), &LINK(dir), &tmat);
      adjoint_f(&tmat, &tmat2);
      sum_mat_f(&tmat, &FS(component));
      dif_mat_f(&tmat2, &FS(component));
    }
    cleanup_gather(mtag);

    // Make traceless
    FORALLSITES(i, s) {
      tc = trace_f(&FS(component));
      CMULREAL(tc, one_ov_N, tc);
      for (j = 0; j < NCOL; j++)
        CSUB(FS(component).e[j][j], tc, FS(component).e[j][j]);
    }
  }
}
// -----------------------------------------------------------------
