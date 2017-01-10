/************* mswvb_gamma_r.c **************************/
/* 
  Multiply a "Wilson matrix" (spin_wilson_vector) by a gamma matrix
  acting on the column index
  (This is the second index, or equivalently, multiplication on the right)
  usage:  mb_gamma_r(src, dest, dir)
  spin_wilson_vector *src,*dest;
  int dir;    dir = XUP, YUP, ZUP, TUP or GAMMAFIVE
*/
#include <stdio.h>
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/dirs.h"

void mult_swv_by_gamma_r(spin_wilson_vector *src, spin_wilson_vector *dest, 
                         int dir) {

  register int i, s1;
  switch(dir) {
    case XUP:
      for (i=0;i<DIMF;i++)for (s1=0;s1<4;s1++) {
        CMUL_MINUS_I(src->d[s1].d[3].c[i], dest->d[s1].d[0].c[i]);
        CMUL_MINUS_I(src->d[s1].d[2].c[i], dest->d[s1].d[1].c[i]);
        CMUL_I(src->d[s1].d[1].c[i], dest->d[s1].d[2].c[i]);
        CMUL_I(src->d[s1].d[0].c[i], dest->d[s1].d[3].c[i]);
      }
      break;
    case YUP:
      for (i=0;i<DIMF;i++)for (s1=0;s1<4;s1++) {
        CNEGATE(src->d[s1].d[3].c[i], dest->d[s1].d[0].c[i]);
        CCOPY(src->d[s1].d[2].c[i], dest->d[s1].d[1].c[i]);
        CCOPY(src->d[s1].d[1].c[i], dest->d[s1].d[2].c[i]);
        CNEGATE(src->d[s1].d[0].c[i], dest->d[s1].d[3].c[i]);
      }
      break;
    case ZUP:
      for (i=0;i<DIMF;i++)for (s1=0;s1<4;s1++) {
        CMUL_MINUS_I(src->d[s1].d[2].c[i], dest->d[s1].d[0].c[i]);
        CMUL_I(src->d[s1].d[3].c[i], dest->d[s1].d[1].c[i]);
        CMUL_I(src->d[s1].d[0].c[i], dest->d[s1].d[2].c[i]);
        CMUL_MINUS_I(src->d[s1].d[1].c[i], dest->d[s1].d[3].c[i]);
      }
      break;
    case TUP:
      for (i=0;i<DIMF;i++)for (s1=0;s1<4;s1++) {
        CCOPY(src->d[s1].d[2].c[i], dest->d[s1].d[0].c[i]);
        CCOPY(src->d[s1].d[3].c[i], dest->d[s1].d[1].c[i]);
        CCOPY(src->d[s1].d[0].c[i], dest->d[s1].d[2].c[i]);
        CCOPY(src->d[s1].d[1].c[i], dest->d[s1].d[3].c[i]);
      }
      break;
    case GAMMAFIVE:
      for (i=0;i<DIMF;i++)for (s1=0;s1<4;s1++) {
        CCOPY(src->d[s1].d[0].c[i], dest->d[s1].d[0].c[i]);
        CCOPY(src->d[s1].d[1].c[i], dest->d[s1].d[1].c[i]);
        CNEGATE(src->d[s1].d[2].c[i], dest->d[s1].d[2].c[i]);
        CNEGATE(src->d[s1].d[3].c[i], dest->d[s1].d[3].c[i]);
      }
      break;
    default:
      printf("BAD CALL TO MULT_BY_GAMMA_RIGHT()\n");
  }
}
