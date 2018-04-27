// -----------------------------------------------------------------
// Blocking function for MCRG-blocked measurements
// Handles "neighboring" sites separated by 2^block links
// Use tempmatf and tempmatf2 for temporary storage
#include "wflow_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Put result into tempmatf2
void staple_mcrg(int dir, int block) {
  register int i, dir2;
  register site *s;
  int j, bl, start, disp[4];    // Displacement vector for general gather
  msg_tag *tag0, *tag1, *tag2, *tag3, *tag4, *tag5, *tag6;
  matrix_f tmat, tmat2;

  bl = 1;
  for (j = 1; j < block; j++)
    bl *= 2;                    // Block size

  start = 1;                    // Indicates staple sum not initialized
  // Loop over other directions
  FORALLUPDIR(dir2) {
    if (dir2 == dir)
      continue;

    // Get linkf[dir2] from direction 2 * dir
    clear_disp(disp);
    disp[dir] = 2 * bl;
    tag0 = start_general_gather_site(F_OFFSET(linkf[dir2]),
                                     sizeof(matrix_f), disp,
                                     EVENANDODD, gen_pt[0]);
    wait_general_gather(tag0);

    // Get linkf[dir] from direction dir2
    clear_disp(disp);
    disp[dir2] = bl;
    tag1 = start_general_gather_site(F_OFFSET(linkf[dir]),
                                     sizeof(matrix_f), disp,
                                     EVENANDODD, gen_pt[1]);
    wait_general_gather(tag1);

    // Get linkf[dir] from direction dir + dir2
    clear_disp(disp);
    disp[dir] = bl;
    disp[dir2] = bl;
    tag2 = start_general_gather_site(F_OFFSET(linkf[dir]),
                                     sizeof(matrix_f), disp,
                                     EVENANDODD, gen_pt[2]);
    wait_general_gather(tag2);

    // Get linkf[dir2] from direction -dir2
    clear_disp(disp);
    disp[dir2] = -bl;
    tag3 = start_general_gather_site(F_OFFSET(linkf[dir2]),
                                     sizeof(matrix_f), disp,
                                     EVENANDODD, gen_pt[3]);
    wait_general_gather(tag3);

    // Get linkf[dir] from direction -dir2
    clear_disp(disp);
    disp[dir2] = -bl;
    tag4 = start_general_gather_site(F_OFFSET(linkf[dir]),
                                     sizeof(matrix_f), disp,
                                     EVENANDODD, gen_pt[4]);
    wait_general_gather(tag4);

    // Get linkf[dir] from direction dir - dir2
    clear_disp(disp);
    disp[dir] = bl;
    disp[dir2]= -bl;
    tag5 = start_general_gather_site(F_OFFSET(linkf[dir]),
                                     sizeof(matrix_f), disp,
                                     EVENANDODD, gen_pt[5]);
    wait_general_gather(tag5);

    // Get linkf[dir2] from direction 2 * dir - dir2
    clear_disp(disp);
    disp[dir] = 2 * bl;
    disp[dir2] = -bl;

    tag6 = start_general_gather_site(F_OFFSET(linkf[dir2]),
                                     sizeof(matrix_f), disp,
                                     EVENANDODD, gen_pt[6]);
    wait_general_gather(tag6);

    // Upper staple
    if (start) {          // The first contribution to the staple
      FORALLSITES(i, s) {
        mult_nn_f(&(s->linkf[dir2]), (matrix_f *)(gen_pt[1][i]), &tmat);
        mult_nn_f(&tmat, (matrix_f *)(gen_pt[2][i]), &tmat2);
        mult_na_f(&tmat2, (matrix_f *)(gen_pt[0][i]), &(tempmatf2[i]));
      }
      start = 0;
    }
    else {
      FORALLSITES(i, s) {
        mult_nn_f(&(s->linkf[dir2]), (matrix_f *)(gen_pt[1][i]), &tmat);
        mult_nn_f(&tmat, (matrix_f *)(gen_pt[2][i]), &tmat2);
        mult_na_sum_f(&tmat2, (matrix_f *)(gen_pt[0][i]), &(tempmatf2[i]));
      }
    }
    cleanup_general_gather(tag0);
    cleanup_general_gather(tag1);
    cleanup_general_gather(tag2);

    // Lower staple
    FORALLSITES(i, s) {
      mult_an_f((matrix_f *)(gen_pt[3][i]), (matrix_f *)(gen_pt[4][i]), &tmat);
      mult_nn_f(&tmat, (matrix_f *)(gen_pt[5][i]), &tmat2);
      mult_nn_sum_f(&tmat2, (matrix_f *)(gen_pt[6][i]), &(tempmatf2[i]));
    }
    cleanup_general_gather(tag3);
    cleanup_general_gather(tag4);
    cleanup_general_gather(tag5);
    cleanup_general_gather(tag6);
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Closely resembles ../generic_nhyp/block_nhyp.c:block_nhyp3()
// If num == 0, do only center link, not the full staple
void block_mcrg(int num, int block) {
  register int dir, i;
  register site *s;
  int disp[4], j;
  Real f[3], tr, tr2;
  complex tc;
  msg_tag *tag0;
  matrix_f tmat, Omega, eQ, Id, Q, Q2;

  FORALLUPDIR(dir) {
    // First the central spine of the blocking
    clear_disp(disp);
    disp[dir] = 1;
    for (j = 1; j < block; j++)
      disp[dir] *= 2;             // Block size

    tag0 = start_general_gather_site(F_OFFSET(linkf[dir]),
                                     sizeof(matrix_f), disp,
                                     EVENANDODD, gen_pt[0]);
    wait_general_gather(tag0);
    FORALLSITES(i, s)
      mult_nn_f(&(s->linkf[dir]), (matrix_f *)(gen_pt[0][i]), &(tempmatf[i]));

    if (num == 0) {               // Do only center link
      FORALLSITES(i, s)
        mat_copy_f(&(tempmatf[i]), &(s->linkf[dir + 4]));
    }
    else {                        // Do the full staple
      staple_mcrg(dir, block);    // Puts result into tempmatf2

      tr = 1.0 - alpha_smear[num];
      tr2 = alpha_smear[num] / 6.0;
      FORALLSITES(i, s) {
        // Make Omega
        scalar_mult_mat_f(&(tempmatf[i]), tr, &Q);
        scalar_mult_add_mat_f(&Q, &(tempmatf2[i]), tr2, &Omega);
        mult_an_f(&Omega, &Omega, &Q);

        // IR stabilization regulator set in defines.h
        scalar_add_diag_f(&Q, IR_STAB);
#ifndef NHYP_DEBUG
        compute_fhb(&Q, f, NULL, 0);
#else
        compute_fhb(&Omega, &Q, f, NULL, 0);
#endif

        // Compute Q**2
        mult_nn_f(&Q, &Q, &Q2);

        // Compute Q^(-1/2) via Eq. (3.8)
        tc = cmplx(f[0], 0);
        diag(&Id, &tc);
        scalar_mult_add_mat_f(&Id, &Q, f[1], &tmat);
        scalar_mult_add_mat_f(&tmat, &Q2, f[2], &eQ);

        // Multiply Omega by eQ = (Omega^\dagger Omega)^(-1/2)
        mult_nn_f(&Omega, &eQ, &tmat);
        mat_copy_f(&tmat, &(s->linkf[dir + 4]));
      }
    }
  }

  FORALLSITES(i, s) {
    FORALLUPDIR(dir)
      mat_copy_f(&(s->linkf[dir + 4]), &(s->linkf[dir]));
  }
}
// -----------------------------------------------------------------
