// -----------------------------------------------------------------
// Compute sum over spins of wilson_vector outer product
// c_ij <- sum(a_i * bdag_j)
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void projector_w(wilson_vector *a, wilson_vector *b, matrix *c) {
  register int i, j, k;
#ifndef FAST
  for (i = 0; i < DIMF; i++) {
    for (j = 0; j < DIMF; j++) {
      CMUL_J(a->d[0].c[i], b->d[0].c[j], c->e[i][j]);
      for (k = 1; k < 4; k++)
        CMUL_JSUM(a->d[k].c[i], b->d[k].c[j], c->e[i][j]);
    }
  }

#else  // FAST version for NCOL = DIMF = 3
  register Real tmp_r, tmp_i, tmp2;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      tmp_r = a->d[0].c[i].real * b->d[0].c[j].real;
      tmp2 = a->d[0].c[i].imag * b->d[0].c[j].imag; tmp_r = tmp_r + tmp2;
      tmp_i = a->d[0].c[i].imag * b->d[0].c[j].real;
      tmp2 = a->d[0].c[i].real * b->d[0].c[j].imag; tmp_i = tmp_i - tmp2;
      for (k = 1; k < 4; k++) {
        tmp2 = a->d[k].c[i].real * b->d[k].c[j].real; tmp_r = tmp_r + tmp2;
        tmp2 = a->d[k].c[i].imag * b->d[k].c[j].imag; tmp_r = tmp_r + tmp2;
        tmp2 = a->d[k].c[i].imag * b->d[k].c[j].real; tmp_i = tmp_i + tmp2;
        tmp2 = a->d[k].c[i].real * b->d[k].c[j].imag; tmp_i = tmp_i - tmp2;
      }
      c->e[i][j].real = tmp_r;
      c->e[i][j].imag = tmp_i;
    }
  }
#endif
}
// -----------------------------------------------------------------
