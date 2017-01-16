// Evaluate the Polyakov loop using general_gathers
// Assume that nt is even
// Stripped PLOOPDIST, which was wrong
/* Reports the distribution of ploop values as a function of each x y z */
/* Similar to ploop3.c, except doesn't use tmat */
// It does use tempmatf and tempmatf2
#include "generic_includes.h"

complex ploop() {
  register int i, t;
  register site *s;
  int d[4] = {0, 0, 0, 0};
  complex sum = cmplx(0.0, 0.0), plp;
  msg_tag *tag;

  /* First multiply the link on every even site by the link above it */
  /* We will compute the Polyakov loop "at" the even sites in the
     first two time slices. */
  tag = start_gather_site(F_OFFSET(linkf[TUP]), sizeof(su3_matrix_f),
                          TUP, EVEN, gen_pt[0]);
  wait_gather(tag);
  FOREVENSITES(i, s) {
    mult_su3_nn_f(&(s->linkf[TUP]), (su3_matrix_f *)gen_pt[0][i],
                  &(tempmatf[i]));
  }
  cleanup_gather(tag);

  for(t = 2; t < nt; t += 2) {
    d[TUP] = t; /* distance from which to gather */
    tag = start_general_gather_field(tempmatf, sizeof(su3_matrix_f),
                                     d, EVEN, gen_pt[0]);
    wait_general_gather(tag);
    FOREVENSITES(i, s) {
      if (s->t > 1)
        continue;  /* only compute on first two slices */
      mult_su3_nn_f(&(tempmatf[i]), (su3_matrix_f *)gen_pt[0][i],
                    &(tempmatf2[i]));
      su3mat_copy_f(&(tempmatf2[i]), &(tempmatf[i]));
      /* We overwrite tempmat1 on the first two time slices,
         leaving the others undisturbed so we can still gather
         them. */
    }
    cleanup_general_gather(tag);
  }
  FOREVENSITES(i, s) {
    if (s->t > 1)
      continue;
    plp = trace_su3_f(&(tempmatf[i]));
    CSUM(sum, plp);
  }
  g_complexsum(&sum);
  plp.real = sum.real / ((Real)(nx * ny * nz));
  plp.imag = sum.imag / ((Real)(nx * ny * nz));
  return plp;
}
// -----------------------------------------------------------------
