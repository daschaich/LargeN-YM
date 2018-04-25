// -----------------------------------------------------------------
// Translate two-index-antisymmetric force term to fundamental rep
// Closely follows make_fermion_rep_matrix
// Basis states are
//   [01 02 03 ... 12 13 ... 23 ...]
#include "cl_dyn_includes.h"

#if DIMF != NCOL * (NCOL - 1) / 2
  #error "Wrong version of chain_rule!"
#endif
#if FREP != ANTISYMMETRIC2
  #error "Wrong version of chain_rule!"
#endif

void chain_rule(matrix_f *sigmaf, matrix *sigma, matrix_f *gaugelinkf) {
  int i, j, k, l;
  int ij = 0, kl;       // Compound indices of the fermion rep
  complex y;

  clear_mat_f(sigmaf);
  for (i = 0; i < NCOL; i++) {
    for (j = i + 1; j < NCOL; j++) {                      // j > i
      kl = 0;
      for (k = 0; k < NCOL; k++) {
        for (l = k + 1; l < NCOL; l++) {                  // l > k
          CMUL(gaugelinkf->e[i][k], sigma->e[kl][ij], y);
          CSUM(sigmaf->e[l][j], y);
          CMUL(gaugelinkf->e[j][k], sigma->e[kl][ij], y);
          CSUB(sigmaf->e[l][i], y, sigmaf->e[l][i]);
          CMUL(gaugelinkf->e[i][l], sigma->e[kl][ij], y);
          CSUB(sigmaf->e[k][j], y, sigmaf->e[k][j]);
          CMUL(gaugelinkf->e[j][l], sigma->e[kl][ij], y);
          CSUM(sigmaf->e[k][i], y);
          kl++;
        }
      }
      ij++;
    }
  }
}
// -----------------------------------------------------------------
