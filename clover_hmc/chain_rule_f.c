// -----------------------------------------------------------------
// The fundamental chain rule is trivial
// Included for consistency with prototype
#include "cl_dyn_includes.h"

#if DIMF != NCOL
  #error "Wrong version of chain_rule!"
#endif
#if FREP != fundamental
  #error "Wrong version of chain_rule!"
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void chain_rule(su3_matrix_f *sigmaf, su3_matrix *sigma,
                su3_matrix_f *gaugelinkf) {

  int i, j;

  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++)
      sigmaf->e[i][j] = sigma->e[i][j];
  }
}
// -----------------------------------------------------------------
