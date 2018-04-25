// -----------------------------------------------------------------
// Translate two-index-symmetric force term to fundamental rep
// Closely follows make_fermion_rep_matrix
// Basis states are
//   [00 01 02 03 ... 11 12 13 ... 22 23 ...]
#include "cl_dyn_includes.h"

#if DIMF != NCOL * (NCOL + 1) / 2
  #error "Wrong version of chain_rule!"
#endif
#if FREP != SYMMETRIC2
  #error "Wrong version of chain_rule!"
#endif

void chain_rule(matrix_f *sigmaf, matrix *sigma, matrix_f *gaugelinkf) {
  int i, j, k, l;
  int ij = 0, kl;       // Compound indices of the fermion rep
  complex x, y;
  double root2 = sqrt(2.0);

  clear_mat_f(sigmaf);
  for (i = 0; i < NCOL; i++) {
    kl = 0;                               // j = i
    for (k = 0; k < NCOL; k++) {
      CMUL(gaugelinkf->e[i][k], sigma->e[kl][ij], x);     // l = k
      CMULREAL(x, 2.0, y);
      CSUM(sigmaf->e[k][i], y);
      kl++;
      for (l = k + 1; l < NCOL; l++) {                    // l > k
        CMUL(gaugelinkf->e[i][k], sigma->e[kl][ij], x);
        CMULREAL(x, root2, y);
        CSUM(sigmaf->e[l][i], y);
        CMUL(gaugelinkf->e[i][l], sigma->e[kl][ij], x);
        CMULREAL(x, root2, y);
        CSUM(sigmaf->e[k][i], y);
        kl++;
      }
    }

    ij++;
    for (j = i + 1; j < NCOL; j++) {      // j > i
      kl = 0;
      for (k = 0; k < NCOL; k++) {
        CMUL(gaugelinkf->e[i][k], sigma->e[kl][ij], x);   // l = k
        CMULREAL(x, root2, y);
        CSUM(sigmaf->e[k][j], y);
        CMUL(gaugelinkf->e[j][k], sigma->e[kl][ij], x);
        CMULREAL(x, root2, y);
        CSUM(sigmaf->e[k][i], y);
        kl++;
        for (l = k + 1; l < NCOL; l++) {                  // l > k
          CMUL(gaugelinkf->e[i][k], sigma->e[kl][ij], y);
          CSUM(sigmaf->e[l][j], y);
          CMUL(gaugelinkf->e[j][k], sigma->e[kl][ij], y);
          CSUM(sigmaf->e[l][i], y);
          CMUL(gaugelinkf->e[i][l], sigma->e[kl][ij], y);
          CSUM(sigmaf->e[k][j], y);
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
