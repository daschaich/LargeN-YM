// -----------------------------------------------------------------
// Compute and compress the traceless anti-hermitian part of a matrix
// dest = 0.5 * (src - src^dag) - Tr[0.5 * (src - src^dag)]
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void make_anti_hermitian(matrix_f *src, anti_hermitmat *dest) {
  Real tr;

#ifndef FAST
  int i, j, index = 0;

  tr = src->e[0][0].imag;
  for (i = 1; i < NCOL; i++)
    tr += src->e[i][i].imag;
  tr /= (Real)NCOL;

  for (i = 0; i < NCOL; i++)
    dest->im_diag[i] = src->e[i][i].imag - tr;

  for (i = 0; i < NCOL; i++) {
    for (j = i + 1; j < NCOL; j++) {
      dest->m[index].real = 0.5 * (src->e[i][j].real - src->e[j][i].real);
      dest->m[index].imag = 0.5 * (src->e[i][j].imag + src->e[j][i].imag);
      index++;
    }
  }

#else  // FAST version for NCOL = DIMF = 3
  Real tr2;
  tr = src->e[0][0].imag + src->e[1][1].imag;
  tr2 = tr + src->e[2][2].imag;
  tr = tr2 * 0.333333333333333333;
  dest->im_diag[0] = src->e[0][0].imag - tr;
  dest->im_diag[1] = src->e[1][1].imag - tr;
  dest->im_diag[2] = src->e[2][2].imag - tr;

  tr = src->e[0][1].real - src->e[1][0].real;
  dest->m[0].real = tr * 0.5;
  tr = src->e[0][1].imag + src->e[1][0].imag;
  dest->m[0].imag = tr * 0.5;

  tr = src->e[0][2].real - src->e[2][0].real;
  dest->m[1].real = tr * 0.5;
  tr = src->e[0][2].imag + src->e[2][0].imag;
  dest->m[1].imag = tr * 0.5;

  tr = src->e[1][2].real - src->e[2][1].real;
  dest->m[2].real = tr * 0.5;
  tr = src->e[1][2].imag + src->e[2][1].imag;
  dest->m[2].imag = tr * 0.5;
#endif
}
// -----------------------------------------------------------------
