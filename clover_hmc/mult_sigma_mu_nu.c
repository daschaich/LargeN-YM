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

  register int i;
  switch(mu) {
    case XUP:
      switch(nu) {
        case YUP:
          for (i = 0; i < DIMF; i++) {
            CMUL_MINUS_I(src->d[0].c[i], dest->d[0].c[i]);
            CMUL_I(src->d[1].c[i], dest->d[1].c[i]);
            CMUL_MINUS_I(src->d[2].c[i], dest->d[2].c[i]);
            CMUL_I(src->d[3].c[i], dest->d[3].c[i]);
          }
          break;
        case ZUP:
          for (i = 0; i < DIMF; i++) {
            CNEGATE(src->d[1].c[i], dest->d[0].c[i]);
            CCOPY(src->d[0].c[i], dest->d[1].c[i]);
            CNEGATE(src->d[3].c[i], dest->d[2].c[i]);
            CCOPY(src->d[2].c[i], dest->d[3].c[i]);
          }
          break;
        case TUP:
          for (i = 0; i < DIMF; i++) {
            CMUL_I(src->d[1].c[i], dest->d[0].c[i]);
            CMUL_I(src->d[0].c[i], dest->d[1].c[i]);
            CMUL_MINUS_I(src->d[3].c[i], dest->d[2].c[i]);
            CMUL_MINUS_I(src->d[2].c[i], dest->d[3].c[i]);
          }
          break;
        default:
          printf("BAD CALL BY mult_sigma_mu_nu\n");
      }
      break;
    case YUP:
      switch(nu) {
        case XUP:
          for (i = 0; i < DIMF; i++) {
            CMUL_I(src->d[0].c[i], dest->d[0].c[i]);
            CMUL_MINUS_I(src->d[1].c[i], dest->d[1].c[i]);
            CMUL_I(src->d[2].c[i], dest->d[2].c[i]);
            CMUL_MINUS_I(src->d[3].c[i], dest->d[3].c[i]);
          }
          break;
        case ZUP:
          for (i = 0; i < DIMF; i++) {
            CMUL_MINUS_I(src->d[1].c[i], dest->d[0].c[i]);
            CMUL_MINUS_I(src->d[0].c[i], dest->d[1].c[i]);
            CMUL_MINUS_I(src->d[3].c[i], dest->d[2].c[i]);
            CMUL_MINUS_I(src->d[2].c[i], dest->d[3].c[i]);
          }
          break;
        case TUP:
          for (i = 0; i < DIMF; i++) {
            CNEGATE(src->d[1].c[i], dest->d[0].c[i]);
            CCOPY(src->d[0].c[i], dest->d[1].c[i]);
            CCOPY(src->d[3].c[i], dest->d[2].c[i]);
            CNEGATE(src->d[2].c[i], dest->d[3].c[i]);
          }
          break;
        default:
          printf("BAD CALL BY mult_sigma_mu_nu\n");
      }
      break;
    case ZUP:
      switch(nu) {
        case XUP:
          for (i = 0; i < DIMF; i++) {
            CCOPY(src->d[1].c[i], dest->d[0].c[i]);
            CNEGATE(src->d[0].c[i], dest->d[1].c[i]);
            CCOPY(src->d[3].c[i], dest->d[2].c[i]);
            CNEGATE(src->d[2].c[i], dest->d[3].c[i]);
          }
          break;
        case YUP:
          for (i = 0; i < DIMF; i++) {
            CMUL_I(src->d[1].c[i], dest->d[0].c[i]);
            CMUL_I(src->d[0].c[i], dest->d[1].c[i]);
            CMUL_I(src->d[3].c[i], dest->d[2].c[i]);
            CMUL_I(src->d[2].c[i], dest->d[3].c[i]);
          }
          break;
        case TUP:
          for (i = 0; i < DIMF; i++) {
            CMUL_I(src->d[0].c[i], dest->d[0].c[i]);
            CMUL_MINUS_I(src->d[1].c[i], dest->d[1].c[i]);
            CMUL_MINUS_I(src->d[2].c[i], dest->d[2].c[i]);
            CMUL_I(src->d[3].c[i], dest->d[3].c[i]);
          }
          break;
        default:
          printf("BAD CALL BY mult_sigma_mu_nu\n");
      }
      break;
    case TUP:
      switch(nu) {
        case XUP:
          for (i = 0; i < DIMF; i++) {
            CMUL_MINUS_I(src->d[1].c[i], dest->d[0].c[i]);
            CMUL_MINUS_I(src->d[0].c[i], dest->d[1].c[i]);
            CMUL_I(src->d[3].c[i], dest->d[2].c[i]);
            CMUL_I(src->d[2].c[i], dest->d[3].c[i]);
          }
          break;
        case YUP:
          for (i = 0; i < DIMF; i++) {
            CCOPY(src->d[1].c[i], dest->d[0].c[i]);
            CNEGATE(src->d[0].c[i], dest->d[1].c[i]);
            CNEGATE(src->d[3].c[i], dest->d[2].c[i]);
            CCOPY(src->d[2].c[i], dest->d[3].c[i]);
          }
          break;
        case ZUP:
          for (i = 0; i < DIMF; i++) {
            CMUL_MINUS_I(src->d[0].c[i], dest->d[0].c[i]);
            CMUL_I(src->d[1].c[i], dest->d[1].c[i]);
            CMUL_I(src->d[2].c[i], dest->d[2].c[i]);
            CMUL_MINUS_I(src->d[3].c[i], dest->d[3].c[i]);
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
