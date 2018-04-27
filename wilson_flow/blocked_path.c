// -----------------------------------------------------------------
// An arbitrary path walker using ordinary gathers
#include "wflow_includes.h"

// Put result in tempmatf2
// Use tempmatf for temporary storage
void blocked_path(int block, int *dir, int *sign, int length) {
  register int i;
  register site *s;
  int j, k, bl;
  msg_tag *mtag;

  bl = 1;
  for (j = 1; j < block; j++)
    bl *= 2;                    // Block size

  if (sign[0] > 0) {
    FORALLSITES(i, s)
      mat_copy_f(&(s->linkf[dir[0]]), &(tempmatf2[i]));

    // Shift block times
    for (k = 0; k < bl; k++) {
      mtag = start_gather_field(tempmatf2, sizeof(matrix_f),
                                OPP_DIR(dir[0]), EVENANDODD, gen_pt[0]);
      wait_gather(mtag);
      FORALLSITES(i, s)
        mat_copy_f((matrix_f *)(gen_pt[0][i]), &(tempmatf[i]));

      cleanup_gather(mtag);
      mtag = start_gather_field(tempmatf, sizeof(matrix_f),
                                OPP_DIR(dir[0]), EVENANDODD, gen_pt[0]);
      wait_gather(mtag);
      FORALLSITES(i, s)
        mat_copy_f((matrix_f *)(gen_pt[0][i]), &(tempmatf2[i]));

      cleanup_gather(mtag);
    }
  }

  if (sign[0] < 0) {
    FORALLSITES(i, s)
      adjoint_f(&(s->linkf[dir[0]]), &(tempmatf2[i]));
  }

  for (j = 1; j < length; j++) {
    if (sign[j] > 0) {
      FORALLSITES(i, s)
        mult_nn_f(&(tempmatf2[i]), &(s->linkf[dir[j]]), &(tempmatf[i]));

      // Shift block times
      for (k = 0; k < bl; k++) {
        mtag = start_gather_field(tempmatf, sizeof(matrix_f),
                                  OPP_DIR(dir[j]), EVENANDODD, gen_pt[0]);
        wait_gather(mtag);
        FORALLSITES(i, s)
          mat_copy_f((matrix_f *)(gen_pt[0][i]), &(tempmatf2[i]));

        cleanup_gather(mtag);
        mtag = start_gather_field(tempmatf2, sizeof(matrix_f),
                                  OPP_DIR(dir[j]), EVENANDODD, gen_pt[0]);
        wait_gather(mtag);
        FORALLSITES(i, s)
          mat_copy_f((matrix_f *)(gen_pt[0][i]), &(tempmatf[i]));

        cleanup_gather(mtag);
      }
      FORALLSITES(i, s)
        mat_copy_f(&(tempmatf[i]), &(tempmatf2[i]));
    }

    if (sign[j] < 0) {
      // Shift block times
      for (k = 0; k < bl; k++) {
        mtag = start_gather_field(tempmatf2, sizeof(matrix_f),
                                  dir[j], EVENANDODD, gen_pt[0]);
        wait_gather(mtag);
        FORALLSITES(i, s)
          mat_copy_f((matrix_f *)(gen_pt[0][i]), &(tempmatf[i]));

        cleanup_gather(mtag);
        mtag = start_gather_field(tempmatf, sizeof(matrix_f),
                                  dir[j], EVENANDODD, gen_pt[0]);
        wait_gather(mtag);
        FORALLSITES(i, s)
          mat_copy_f((matrix_f *)(gen_pt[0][i]), &(tempmatf2[i]));

        cleanup_gather(mtag);
      }
      FORALLSITES(i, s)
        mat_copy_f(&(tempmatf2[i]), &(tempmatf[i]));

      FORALLSITES(i, s)
        mult_na_f(&(tempmatf[i]), &(s->linkf[dir[j]]), &(tempmatf2[i]));
    }
  }
}
// -----------------------------------------------------------------
