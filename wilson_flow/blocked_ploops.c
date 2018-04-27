// -----------------------------------------------------------------
// Evaluate the Polyakov loop or Wilson line in the given direction
// Use general_gathers; lattice must be divisible by 2^block in all dirs
// Input "block" reports how many times lattice has been blocked
// Use tempmatf for temporary storage
#include "wflow_includes.h"

complex blocked_ploop(int block, int dir) {
  register int i, k;
  register site *s;
  int bl = 2, d[4], N = nt;
  complex plp = cmplx(0.0, 0.0);
  msg_tag *tag;

  // Sanity check: reproduce ploop(dir) with this routine
  if (block <= 0)
    bl = 1;

  // Set number of links to stride, bl = 2^block
  // This usage may differ from blocked_path.c, where bl = 2^{block - 1}
  for (k = 1; k < block; k++)
    bl *= 2;

  FORALLUPDIR(i)
    d[i] = 0;

  switch(dir) {
    case XUP: N = nx; break;
    case YUP: N = ny; break;
    case ZUP: N = nz; break;
    case TUP: N = nt; break;
    default:
      node0_printf("ERROR: Unrecognized direction in blocked_ploop\n");
      terminate(1);
  }

  // Copy links to tempmatf
  FORALLSITES(i, s)
    mat_copy_f(&(s->linkf[dir]), &(tempmatf[i]));

  // Compute the bl-strided Polyakov loop "at" ALL the sites
  // on the first 2^block timeslices
  for (k = bl; k < N; k += bl) {
    d[dir] = k;                     // Distance from which to gather
    tag = start_general_gather_field(tempmatf, sizeof(matrix_f),
                                     d, EVENANDODD, gen_pt[0]);
    wait_general_gather(tag);
    FORALLSITES(i, s) {
      // Overwrite tempmatf on the first bl slices
      // Leave other links undisturbed so we can still gather them
      switch(dir) {
        case XUP: if (s->x >= bl) continue; break;
        case YUP: if (s->y >= bl) continue; break;
        case ZUP: if (s->z >= bl) continue; break;
        case TUP: if (s->t >= bl) continue; break;
      }
      mult_nn_f(&(tempmatf[i]), (matrix_f *)gen_pt[0][i], &(tempmatf2[i]));
      mat_copy_f(&(tempmatf2[i]), &(tempmatf[i]));
    }
    cleanup_general_gather(tag);
  }
  FORALLSITES(i, s) {
    switch(dir) {
      case XUP: if (s->x >= bl) continue; break;
      case YUP: if (s->y >= bl) continue; break;
      case ZUP: if (s->z >= bl) continue; break;
      case TUP: if (s->t >= bl) continue; break;
    }
    trace_sum_f(&(tempmatf[i]), &plp);     // Running complex sum
  }

  // Average all the loops we just calculated
  g_complexsum(&plp);
  plp.real *= N / ((double)(volume * bl));
  plp.imag *= N / ((double)(volume * bl));
  return plp;
}
// -----------------------------------------------------------------
