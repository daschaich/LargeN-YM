// An arbitrary path walker for fundamental links
// Uses ordinary gathers, accumulates in tempmatf
#include "cl_dyn_includes.h"

void path(int *dir, int *sign, int length) {
  register int i, j;
  register site *s;
  msg_tag *mtag0;

  /* j=0 */
  if (sign[0] > 0) {
    mtag0 = start_gather_site(F_OFFSET(linkf[dir[0]]), sizeof(su3_matrix_f),
                              OPP_DIR(dir[0]), EVENANDODD, gen_pt[0]);
    wait_gather(mtag0);

    FORALLSITES(i, s)
      su3mat_copy_f((su3_matrix_f *)(gen_pt[0][i]), &(tempmatf[i]));

    cleanup_gather(mtag0);
  }

  if (sign[0] < 0) {
    FORALLSITES(i, s)
      su3_adjoint_f(&(s->linkf[dir[0]]), &(tempmatf[i]));
  }


  for (j = 1; j < length; j++) {
    if (sign[j] > 0) {

      FORALLSITES(i, s)
        mult_su3_nn_f(&(tempmatf[i]), &(s->linkf[dir[j]]), &(tempmatf2[i]));

      mtag0 = start_gather_field(tempmatf2, sizeof(su3_matrix_f),
                                 OPP_DIR(dir[j]), EVENANDODD, gen_pt[0]);

      wait_gather(mtag0);
      FORALLSITES(i, s)
        su3mat_copy_f((su3_matrix_f *)(gen_pt[0][i]), &(tempmatf[i]));

      cleanup_gather(mtag0);
    }

    if (sign[j] < 0) {
      mtag0 = start_gather_field(tempmatf, sizeof(su3_matrix_f),
                                 dir[j], EVENANDODD, gen_pt[1]);

      wait_gather(mtag0);
      FORALLSITES(i, s) {
        mult_su3_na_f((su3_matrix_f *)(gen_pt[1][i]), &(s->linkf[dir[j]]),
                      &(tempmatf2[i]));
      }

      FORALLSITES(i, s)
        su3mat_copy_f(&(tempmatf2[i]), &(tempmatf[i]));

      cleanup_gather(mtag0);
    }
  }
}
