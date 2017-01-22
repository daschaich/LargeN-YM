/************* mb_gamma_r.c  (in su3.a) **************************/
/*
  Multiply a Wilson matrix by a gamma matrix acting on the column index
  (This is the second index, or equivalently, multiplication on the right)
  usage:   mult_by_gamma_right wilson_matrix *src,  wilson_matrix *dest,
  int dir)
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

void mult_by_gamma_right(wilson_matrix *src, wilson_matrix *dest, int dir) {
  register int i, c1, s1;

  switch(dir) {
    case XUP:
      for (i = 0; i < DIMF; i++)for (s1=0;s1<4;s1++)for (c1=0;c1<DIMF;c1++) {
        CMUL_MINUS_I(src->d[s1].c[c1].d[3].c[i], dest->d[s1].c[c1].d[0].c[i]);
        CMUL_MINUS_I(src->d[s1].c[c1].d[2].c[i], dest->d[s1].c[c1].d[1].c[i]);
        CMUL_I(src->d[s1].c[c1].d[1].c[i], dest->d[s1].c[c1].d[2].c[i]);
        CMUL_I(src->d[s1].c[c1].d[0].c[i], dest->d[s1].c[c1].d[3].c[i]);
      }
      break;
    case YUP:
      for (i = 0; i < DIMF; i++)for (s1=0;s1<4;s1++)for (c1=0;c1<DIMF;c1++) {
        CNEGATE(src->d[s1].c[c1].d[3].c[i], dest->d[s1].c[c1].d[0].c[i]);
        CCOPY(src->d[s1].c[c1].d[2].c[i], dest->d[s1].c[c1].d[1].c[i]);
        CCOPY(src->d[s1].c[c1].d[1].c[i], dest->d[s1].c[c1].d[2].c[i]);
        CNEGATE(src->d[s1].c[c1].d[0].c[i], dest->d[s1].c[c1].d[3].c[i]);
      }
      break;
    case ZUP:
      for (i = 0; i < DIMF; i++)for (s1=0;s1<4;s1++)for (c1=0;c1<DIMF;c1++) {
        CMUL_MINUS_I(src->d[s1].c[c1].d[2].c[i], dest->d[s1].c[c1].d[0].c[i]);
        CMUL_I(src->d[s1].c[c1].d[3].c[i], dest->d[s1].c[c1].d[1].c[i]);
        CMUL_I(src->d[s1].c[c1].d[0].c[i], dest->d[s1].c[c1].d[2].c[i]);
        CMUL_MINUS_I(src->d[s1].c[c1].d[1].c[i], dest->d[s1].c[c1].d[3].c[i]);
      }
      break;
    case TUP:
      for (i = 0; i < DIMF; i++)for (s1=0;s1<4;s1++)for (c1=0;c1<DIMF;c1++) {
        CCOPY(src->d[s1].c[c1].d[2].c[i], dest->d[s1].c[c1].d[0].c[i]);
        CCOPY(src->d[s1].c[c1].d[3].c[i], dest->d[s1].c[c1].d[1].c[i]);
        CCOPY(src->d[s1].c[c1].d[0].c[i], dest->d[s1].c[c1].d[2].c[i]);
        CCOPY(src->d[s1].c[c1].d[1].c[i], dest->d[s1].c[c1].d[3].c[i]);
      }
      break;
    case GAMMAFIVE:
      for (i = 0; i < DIMF; i++)for (s1=0;s1<4;s1++)for (c1=0;c1<DIMF;c1++) {
        CCOPY(src->d[s1].c[c1].d[0].c[i], dest->d[s1].c[c1].d[0].c[i]);
        CCOPY(src->d[s1].c[c1].d[1].c[i], dest->d[s1].c[c1].d[1].c[i]);
        CNEGATE(src->d[s1].c[c1].d[2].c[i], dest->d[s1].c[c1].d[2].c[i]);
        CNEGATE(src->d[s1].c[c1].d[3].c[i], dest->d[s1].c[c1].d[3].c[i]);
      }
      break;
    default:
      printf("BAD CALL TO MULT_BY_GAMMA_RIGHT()\n");
  }
}

