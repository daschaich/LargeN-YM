// -----------------------------------------------------------------
// The fundamental irrep link field creation is trivial
// Included for consistency with prototype
#include "generic_wilson_includes.h"

#if (DIMF != NCOL)
#error "Wrong version of fermion_rep!"
#endif
#if (FREP != FUNDAMENTAL)
#error "Wrong version of fermion_rep!"
#endif

void make_fermion_rep_matrix(matrix_f *a, matrix *b) {
  int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++)
      b->e[i][j] = a->e[i][j];
  }
}
// -----------------------------------------------------------------
