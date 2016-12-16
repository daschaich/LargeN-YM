/*****************  su3_rdot.c  (in su3.a) ******************************
*									*
* Real su3_rdot( su3_vector *a, su3_vector *b )			*
* return real part of dot product of two su3_vectors			*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

Real su3_rdot( su3_vector *a, su3_vector *b ){

#ifndef NATIVEDOUBLE
register Real temp1,temp2;
#if (DIMF<=6)
    temp2 = a->c[0].real * b->c[0].real;
    temp1 = a->c[0].imag * b->c[0].imag; temp2 += temp1;
    temp1 = a->c[1].real * b->c[1].real; temp2 += temp1;
    temp1 = a->c[1].imag * b->c[1].imag; temp2 += temp1;
#if (DIMF>2)
    temp1 = a->c[2].real * b->c[2].real; temp2 += temp1;
    temp1 = a->c[2].imag * b->c[2].imag; temp2 += temp1;
#if (DIMF>3)
    temp1 = a->c[3].real * b->c[3].real; temp2 += temp1;
    temp1 = a->c[3].imag * b->c[3].imag; temp2 += temp1;
#if (DIMF>4)
    temp1 = a->c[4].real * b->c[4].real; temp2 += temp1;
    temp1 = a->c[4].imag * b->c[4].imag; temp2 += temp1;
#if (DIMF>5)
    temp1 = a->c[5].real * b->c[5].real; temp2 += temp1;
    temp1 = a->c[5].imag * b->c[5].imag; temp2 += temp1;
#endif
#endif
#endif
#endif
#else /* DIMF>6 */
    int i;
    temp2 = a->c[0].real * b->c[0].real;
    temp1 = a->c[0].imag * b->c[0].imag; temp2 += temp1;
    for (i=1;i<DIMF;i++){
        temp1 = a->c[i].real * b->c[i].real; temp2 += temp1;
        temp1 = a->c[i].imag * b->c[i].imag; temp2 += temp1;
    }
#endif

    return(temp2);

#else /* RS6000 version */

  register double ar,ai,br,bi,ss;

  ar=a->c[0].real;  ai=a->c[0].imag;
  br=b->c[0].real;  bi=b->c[0].imag;
  ss = ar*br + ai*bi;

  ar=a->c[1].real;  ai=a->c[1].imag;
  br=b->c[1].real;  bi=b->c[1].imag;
  ss += ar*br + ai*bi;

  ar=a->c[2].real;  ai=a->c[2].imag;
  br=b->c[2].real;  bi=b->c[2].imag;
  ss += ar*br + ai*bi;

  return(ss);

#endif
}
