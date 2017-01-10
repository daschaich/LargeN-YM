// -----------------------------------------------------------------
// Flip the temporal boundary conditions (BCs)
// by multiplying by -1 on the last time-slice time
// "sign" is PLUS or MINUS and says what the BCs should become
// Print warning  if the BCs are already set to "sign"

// This is normal switching for measurement routines
// It is called when generating link from linkf, in fermion_rep
// current_boundary is declared in lattice.h
#include "generic_wilson_includes.h"

void boundary_flip(int sign) {
  register int i,j,k;
  register site *s;

  if (sign == current_boundary) {
    node0_printf("WARNING: You lost track of the boundary conditions!\n");
    return;
  }

  FORALLSITES(i, s) {
    if (s->t != nt - 1)
      continue;  /* hit only last time slice */
    for (j = 0; j < DIMF; j++) {
      for (k = 0; k < DIMF; k++) {
        s->link[TUP].e[j][k].real *= -1.0;
        s->link[TUP].e[j][k].imag *= -1.0;
      }
    }
  }
  current_boundary = sign;
}
// -----------------------------------------------------------------
