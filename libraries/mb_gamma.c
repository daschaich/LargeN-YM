/************* mb_gamma.c  (in su3.a) **************************/
/* 
  Multiply a Wilson vector by a gamma matrix
  usage:  mult_by_gamma(wilson_vector *src, wilson_vector *dest, int dir)
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

void mult_by_gamma(wilson_vector *src, wilson_vector *dest, int dir) {
  register int i;

  switch(dir) {
    case XUP:
      for (i = 0; i < DIMF; i++) {
        CMUL_I(src->d[3].c[i], dest->d[0].c[i]);
        CMUL_I(src->d[2].c[i], dest->d[1].c[i]);
        CMUL_MINUS_I(src->d[1].c[i], dest->d[2].c[i]);
        CMUL_MINUS_I(src->d[0].c[i], dest->d[3].c[i]);
      }
      break;
    case YUP:
      for (i = 0; i < DIMF; i++) {
        CNEGATE(src->d[3].c[i], dest->d[0].c[i]);
        CCOPY(src->d[2].c[i], dest->d[1].c[i]);
        CCOPY(src->d[1].c[i], dest->d[2].c[i]);
        CNEGATE(src->d[0].c[i], dest->d[3].c[i]);
      }
      break;
    case ZUP:
      for (i = 0; i < DIMF; i++) {
        CMUL_I(src->d[2].c[i], dest->d[0].c[i]);
        CMUL_MINUS_I(src->d[3].c[i], dest->d[1].c[i]);
        CMUL_MINUS_I(src->d[0].c[i], dest->d[2].c[i]);
        CMUL_I(src->d[1].c[i], dest->d[3].c[i]);
      }
      break;
    case TUP:
      for (i = 0; i < DIMF; i++) {
        CCOPY(src->d[2].c[i], dest->d[0].c[i]);
        CCOPY(src->d[3].c[i], dest->d[1].c[i]);
        CCOPY(src->d[0].c[i], dest->d[2].c[i]);
        CCOPY(src->d[1].c[i], dest->d[3].c[i]);
      }
      break;
    case GAMMAFIVE:
      for (i = 0; i < DIMF; i++) {
        CCOPY(src->d[0].c[i], dest->d[0].c[i]);
        CCOPY(src->d[1].c[i], dest->d[1].c[i]);
        CNEGATE(src->d[2].c[i], dest->d[2].c[i]);
        CNEGATE(src->d[3].c[i], dest->d[3].c[i]);
      }
      break;
    default:
      printf("BAD CALL TO MULT_BY_GAMMA()\n");
  }
}

