// -----------------------------------------------------------------
// An arbitrary path walker using ordinary gathers
#include "wflow_includes.h"

// Put result in tempmatf2
// Use tempmatf for temporary storage
void path(int *dir, int *sign, int length) {
  register int i;
  register site *s;
  int j;
  msg_tag *mtag;

  if (sign[0] > 0)  {
    mtag = start_gather_site(F_OFFSET(linkf[dir[0]]), sizeof(matrix_f),
                             OPP_DIR(dir[0]), EVENANDODD, gen_pt[0]);
    wait_gather(mtag);
    FORALLSITES(i, s)
      mat_copy_f((matrix_f *)(gen_pt[0][i]), &(tempmatf2[i]));

    cleanup_gather(mtag);
  }

  if (sign[0] < 0) {
    FORALLSITES(i, s)
      adjoint_f(&(s->linkf[dir[0]]), &(tempmatf2[i]));
  }

  for (j = 1; j < length; j++) {
    if (sign[j] > 0) {
      FORALLSITES(i, s)
        mult_nn_f(&(tempmatf2[i]), &(s->linkf[dir[j]]), &(tempmatf[i]));

      mtag = start_gather_field(tempmatf, sizeof(matrix_f),
                                OPP_DIR(dir[j]), EVENANDODD, gen_pt[0]);
      wait_gather(mtag);
      FORALLSITES(i, s)
        mat_copy_f((matrix_f *)(gen_pt[0][i]), &(tempmatf2[i]));

      cleanup_gather(mtag);
    }

    if (sign[j] < 0) {
      mtag = start_gather_field(tempmatf2, sizeof(matrix_f),
                                dir[j], EVENANDODD, gen_pt[1]);
      wait_gather(mtag);
      FORALLSITES(i, s) {
        mult_na_f((matrix_f *)(gen_pt[1][i]), &(s->linkf[dir[j]]),
                  &(tempmatf[i]));
      }

      FORALLSITES(i, s)
        mat_copy_f(&(tempmatf[i]), &(tempmatf2[i]));

      cleanup_gather(mtag);
    }
  }
}
// -----------------------------------------------------------------
