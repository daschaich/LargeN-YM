// -----------------------------------------------------------------
// Create two-index-antisymmetric irrep link field from fundamental link
// Basis states are
//   [01 02 03 ... 12 13 ... 23 ...]
#include "generic_wilson_includes.h"

#if DIMF != NCOL * (NCOL - 1) / 2
  #error "Wrong version of fermion_rep!"
#endif
#if FREP != ANTISYMMETRIC2
  #error "Wrong version of fermion_rep!"
#endif

void make_fermion_rep_matrix(matrix_f *a, matrix *b) {
  int i, j, k, l;
  int ij = 0, kl;       // Compound indices of the fermion rep
  complex x, y;

  for (i = 0; i < NCOL; i++) {
    for (j = i + 1; j <NCOL; j++) {
      kl = 0;
      for (k = 0; k < NCOL; k++) {
        for (l = k + 1; l < NCOL; l++) {
          CMUL(a->e[i][k], a->e[j][l], x);
          CMUL(a->e[i][l], a->e[j][k], y);
          CSUB(x, y, b->e[ij][kl]);
          kl++;
        }
      }
      ij++;
    }
  }
}
// -----------------------------------------------------------------
