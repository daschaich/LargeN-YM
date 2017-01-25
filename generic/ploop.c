// Evaluate the Polyakov loop using general_gathers
// Assume that nt is even
// Stripped PLOOPDIST, which was wrong
// Use tempmatf and tempmatf2 for temporary storage
#include "generic_includes.h"

// Compute the Polyakov loop at the even sites in the first two time-slices
complex ploop() {
  register int i, t;
  register site *s;
  int d[4] = {0, 0, 0, 0};
  complex plp = cmplx(0.0, 0.0);
  msg_tag *tag;

  // First multiply the link on every even site by the link above it
  tag = start_gather_site(F_OFFSET(linkf[TUP]), sizeof(matrix_f),
                          TUP, EVEN, gen_pt[0]);
  wait_gather(tag);
  FOREVENSITES(i, s)
    mult_nn_f(&(s->linkf[TUP]), (matrix_f *)gen_pt[0][i], &(tempmatf[i]));
  cleanup_gather(tag);

  for(t = 2; t < nt; t += 2) {
    d[TUP] = t;     // Distance from which to gather
    tag = start_general_gather_field(tempmatf, sizeof(matrix_f),
                                     d, EVEN, gen_pt[0]);
    wait_general_gather(tag);
    FOREVENSITES(i, s) {
      // Overwrite tempmatf on the first two time slices,
      // leaving the others undisturbed so we can still gather them
      if (s->t > 1)
        continue;
      mult_nn_f(&(tempmatf[i]), (matrix_f *)gen_pt[0][i], &(tempmatf2[i]));
      mat_copy_f(&(tempmatf2[i]), &(tempmatf[i]));
    }
    cleanup_general_gather(tag);
  }
  FOREVENSITES(i, s) {
    if (s->t > 1)
      continue;
    trace_sum_f(&(tempmatf[i]), &plp);
  }
  g_complexsum(&plp);
  plp.real /= (Real)(nx * ny * nz);
  plp.imag /= (Real)(nx * ny * nz);
  return plp;
}
// -----------------------------------------------------------------
