/* Left multiplies a wilson_vector by sigma_mu_nu = gamma_mu gamma_nu where

sigma(XUP,YUP)    sigma(XUP,ZUP)    sigma(XUP,TUP)
  -i  0  0  0    0 -1  0  0    0  i  0  0
   0  i  0  0    1  0  0  0    i  0  0  0
   0  0 -i  0    0  0  0 -1    0  0  0 -i
   0  0  0  i    0  0  1  0    0  0 -i  0

sigma(YUP,ZUP)    sigma(YUP,TUP)    sigma(ZUP,TUP)
   0 -i  0  0    0 -1  0  0    i  0  0  0
  -i  0  0  0    1  0  0  0    0 -i  0  0
   0  0  0 -i    0  0  0  1    0  0 -i  0
   0  0 -i  0    0  0 -1  0    0  0  0  i

and sigma(nu,mu) = -sigma(mu,nu)
*/

#include "cl_dyn_includes.h"

void mult_sigma_mu_nu(wilson_vector *src, wilson_vector *dest,
                      int mu, int nu) {

  register int c2;
  switch(mu) {
    case XUP:
      switch(nu) {
        case YUP:
          for(c2=0;c2<DIMF;c2++) {
            CMUL_MINUS_I(src->d[0].c[c2], dest->d[0].c[c2]);
            CMUL_I(src->d[1].c[c2], dest->d[1].c[c2]);
            CMUL_MINUS_I(src->d[2].c[c2], dest->d[2].c[c2]);
            CMUL_I(src->d[3].c[c2], dest->d[3].c[c2]);
          }
          break;
        case ZUP:
          for(c2=0;c2<DIMF;c2++) {
            CNEGATE(src->d[1].c[c2], dest->d[0].c[c2]);
            CCOPY(src->d[0].c[c2], dest->d[1].c[c2]);
            CNEGATE(src->d[3].c[c2], dest->d[2].c[c2]);
            CCOPY(src->d[2].c[c2], dest->d[3].c[c2]);
          }
          break;
        case TUP:
          for(c2=0;c2<DIMF;c2++) {
            CMUL_I(src->d[1].c[c2], dest->d[0].c[c2]);
            CMUL_I(src->d[0].c[c2], dest->d[1].c[c2]);
            CMUL_MINUS_I(src->d[3].c[c2], dest->d[2].c[c2]);
            CMUL_MINUS_I(src->d[2].c[c2], dest->d[3].c[c2]);
          }
          break;
        default:
          printf("BAD CALL BY mult_sigma_mu_nu\n");
      }
      break;
    case YUP:
      switch(nu) {
        case XUP:
          for(c2=0;c2<DIMF;c2++) {
            CMUL_I(src->d[0].c[c2], dest->d[0].c[c2]);
            CMUL_MINUS_I(src->d[1].c[c2], dest->d[1].c[c2]);
            CMUL_I(src->d[2].c[c2], dest->d[2].c[c2]);
            CMUL_MINUS_I(src->d[3].c[c2], dest->d[3].c[c2]);
          }
          break;
        case ZUP:
          for(c2=0;c2<DIMF;c2++) {
            CMUL_MINUS_I(src->d[1].c[c2], dest->d[0].c[c2]);
            CMUL_MINUS_I(src->d[0].c[c2], dest->d[1].c[c2]);
            CMUL_MINUS_I(src->d[3].c[c2], dest->d[2].c[c2]);
            CMUL_MINUS_I(src->d[2].c[c2], dest->d[3].c[c2]);
          }
          break;
        case TUP:
          for(c2=0;c2<DIMF;c2++) {
            CNEGATE(src->d[1].c[c2], dest->d[0].c[c2]);
            CCOPY(src->d[0].c[c2], dest->d[1].c[c2]);
            CCOPY(src->d[3].c[c2], dest->d[2].c[c2]);
            CNEGATE(src->d[2].c[c2], dest->d[3].c[c2]);
          }
          break;
        default:
          printf("BAD CALL BY mult_sigma_mu_nu\n");
      }
      break;
    case ZUP:
      switch(nu) {
        case XUP:
          for(c2=0;c2<DIMF;c2++) {
            CCOPY(src->d[1].c[c2], dest->d[0].c[c2]);
            CNEGATE(src->d[0].c[c2], dest->d[1].c[c2]);
            CCOPY(src->d[3].c[c2], dest->d[2].c[c2]);
            CNEGATE(src->d[2].c[c2], dest->d[3].c[c2]);
          }
          break;
        case YUP:
          for(c2=0;c2<DIMF;c2++) {
            CMUL_I(src->d[1].c[c2], dest->d[0].c[c2]);
            CMUL_I(src->d[0].c[c2], dest->d[1].c[c2]);
            CMUL_I(src->d[3].c[c2], dest->d[2].c[c2]);
            CMUL_I(src->d[2].c[c2], dest->d[3].c[c2]);
          }
          break;
        case TUP:
          for(c2=0;c2<DIMF;c2++) {
            CMUL_I(src->d[0].c[c2], dest->d[0].c[c2]);
            CMUL_MINUS_I(src->d[1].c[c2], dest->d[1].c[c2]);
            CMUL_MINUS_I(src->d[2].c[c2], dest->d[2].c[c2]);
            CMUL_I(src->d[3].c[c2], dest->d[3].c[c2]);
          }
          break;
        default:
          printf("BAD CALL BY mult_sigma_mu_nu\n");
      }
      break;
    case TUP:
      switch(nu) {
        case XUP:
          for(c2=0;c2<DIMF;c2++) {
            CMUL_MINUS_I(src->d[1].c[c2], dest->d[0].c[c2]);
            CMUL_MINUS_I(src->d[0].c[c2], dest->d[1].c[c2]);
            CMUL_I(src->d[3].c[c2], dest->d[2].c[c2]);
            CMUL_I(src->d[2].c[c2], dest->d[3].c[c2]);
          }
          break;
        case YUP:
          for(c2=0;c2<DIMF;c2++) {
            CCOPY(src->d[1].c[c2], dest->d[0].c[c2]);
            CNEGATE(src->d[0].c[c2], dest->d[1].c[c2]);
            CNEGATE(src->d[3].c[c2], dest->d[2].c[c2]);
            CCOPY(src->d[2].c[c2], dest->d[3].c[c2]);
          }
          break;
        case ZUP:
          for(c2=0;c2<DIMF;c2++) {
            CMUL_MINUS_I(src->d[0].c[c2], dest->d[0].c[c2]);
            CMUL_I(src->d[1].c[c2], dest->d[1].c[c2]);
            CMUL_I(src->d[2].c[c2], dest->d[2].c[c2]);
            CMUL_MINUS_I(src->d[3].c[c2], dest->d[3].c[c2]);
          }
          break;
        default:
          printf("BAD CALL BY mult_sigma_mu_nu\n");
      }
      break;
    default:
      printf("BAD CALL BY mult_sigma_mu_nu\n");
  }

}
