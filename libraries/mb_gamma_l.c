/************* mb_gamma_l.c  (in su3.a) **************************/
/*
  Multiply a Wilson matrix by a gamma matrix acting on the row index
  (This is the first index, or equivalently, multiplication on the left)
  usage: mult_by_gamma_left(wilson_matrix *src,  wilson_matrix *dest, int dir)
  dir = XUP, YUP, ZUP, TUP or GAMMAFIVE

 gamma(XUP) 0  0  0  i
            0  0  i  0
            0 -i  0  0
           -i  0  0  0

 gamma(YUP) 0  0  0 -1
            0  0  1  0
            0  1  0  0
           -1  0  0  0

 gamma(ZUP) 0  0  i  0
            0  0  0 -i
           -i  0  0  0
            0  i  0  0

 gamma(TUP) 0  0  1  0
            0  0  0  1
            1  0  0  0
            0  1  0  0

 gamma(FIVE) 1  0  0  0
             0  1  0  0
             0  0 -1  0
             0  0  0 -1
*/
#include <stdio.h>
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/dirs.h"

void mult_by_gamma_left(wilson_matrix *src, wilson_matrix *dest, int dir) {
  register int i, c2, s2;

  switch(dir) {
    case XUP:
      for (i = 0; i < DIMF; i++)for (s2=0;s2<4;s2++)for (c2=0;c2<DIMF;c2++) {
        CMUL_I(src->d[3].c[i].d[s2].c[c2], dest->d[0].c[i].d[s2].c[c2]);
        CMUL_I(src->d[2].c[i].d[s2].c[c2], dest->d[1].c[i].d[s2].c[c2]);
        CMUL_MINUS_I(src->d[1].c[i].d[s2].c[c2], dest->d[2].c[i].d[s2].c[c2]);
        CMUL_MINUS_I(src->d[0].c[i].d[s2].c[c2], dest->d[3].c[i].d[s2].c[c2]);
      }
      break;
    case YUP:
      for (i = 0; i < DIMF; i++)for (s2=0;s2<4;s2++)for (c2=0;c2<DIMF;c2++) {
        CNEGATE(src->d[3].c[i].d[s2].c[c2], dest->d[0].c[i].d[s2].c[c2]);
        CCOPY(src->d[2].c[i].d[s2].c[c2], dest->d[1].c[i].d[s2].c[c2]);
        CCOPY(src->d[1].c[i].d[s2].c[c2], dest->d[2].c[i].d[s2].c[c2]);
        CNEGATE(src->d[0].c[i].d[s2].c[c2], dest->d[3].c[i].d[s2].c[c2]);
      }
      break;
    case ZUP:
      for (i = 0; i < DIMF; i++)for (s2=0;s2<4;s2++)for (c2=0;c2<DIMF;c2++) {
        CMUL_I(src->d[2].c[i].d[s2].c[c2], dest->d[0].c[i].d[s2].c[c2]);
        CMUL_MINUS_I(src->d[3].c[i].d[s2].c[c2], dest->d[1].c[i].d[s2].c[c2]);
        CMUL_MINUS_I(src->d[0].c[i].d[s2].c[c2], dest->d[2].c[i].d[s2].c[c2]);
        CMUL_I(src->d[1].c[i].d[s2].c[c2], dest->d[3].c[i].d[s2].c[c2]);
      }
      break;
    case TUP:
      for (i = 0; i < DIMF; i++)for (s2=0;s2<4;s2++)for (c2=0;c2<DIMF;c2++) {
        CCOPY(src->d[2].c[i].d[s2].c[c2], dest->d[0].c[i].d[s2].c[c2]);
        CCOPY(src->d[3].c[i].d[s2].c[c2], dest->d[1].c[i].d[s2].c[c2]);
        CCOPY(src->d[0].c[i].d[s2].c[c2], dest->d[2].c[i].d[s2].c[c2]);
        CCOPY(src->d[1].c[i].d[s2].c[c2], dest->d[3].c[i].d[s2].c[c2]);
      }
      break;
    case GAMMAFIVE:
      for (i = 0; i < DIMF; i++)for (s2=0;s2<4;s2++)for (c2=0;c2<DIMF;c2++) {
        CCOPY(src->d[0].c[i].d[s2].c[c2], dest->d[0].c[i].d[s2].c[c2]);
        CCOPY(src->d[1].c[i].d[s2].c[c2], dest->d[1].c[i].d[s2].c[c2]);
        CNEGATE(src->d[2].c[i].d[s2].c[c2], dest->d[2].c[i].d[s2].c[c2]);
        CNEGATE(src->d[3].c[i].d[s2].c[c2], dest->d[3].c[i].d[s2].c[c2]);
      }
      break;
    default:
      printf("BAD CALL TO MULT_BY_GAMMA_LEFT()\n");
  }
}

