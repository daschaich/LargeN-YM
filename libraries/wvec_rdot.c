/*****************  wvec_rdot.c  (in su3.a) ******************************
*									*
* Real wvec_rdot( wilson_vector *a, wilson_vector *b )			*
* return real part of dot product of two wilson_vectors			*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

Real wvec_rdot( wilson_vector *a, wilson_vector *b ){

#ifndef FAST
  register Real temp1,temp2;
  register int i;
#if (DIMF<=6)
  temp2=0.0;
  for(i=0;i<4;i++){
        temp1 = a->d[i].c[0].real * b->d[i].c[0].real; temp2 += temp1;
        temp1 = a->d[i].c[0].imag * b->d[i].c[0].imag; temp2 += temp1;
        temp1 = a->d[i].c[1].real * b->d[i].c[1].real; temp2 += temp1;
        temp1 = a->d[i].c[1].imag * b->d[i].c[1].imag; temp2 += temp1;
#if (DIMF>2)
        temp1 = a->d[i].c[2].real * b->d[i].c[2].real; temp2 += temp1;
        temp1 = a->d[i].c[2].imag * b->d[i].c[2].imag; temp2 += temp1;
#if (DIMF>3)
        temp1 = a->d[i].c[3].real * b->d[i].c[3].real; temp2 += temp1;
        temp1 = a->d[i].c[3].imag * b->d[i].c[3].imag; temp2 += temp1;
#if (DIMF>4)
        temp1 = a->d[i].c[4].real * b->d[i].c[4].real; temp2 += temp1;
        temp1 = a->d[i].c[4].imag * b->d[i].c[4].imag; temp2 += temp1;
#if (DIMF>5)
        temp1 = a->d[i].c[5].real * b->d[i].c[5].real; temp2 += temp1;
        temp1 = a->d[i].c[5].imag * b->d[i].c[5].imag; temp2 += temp1;
#endif
#endif
#endif
#endif
  }
#else /* DIMF>6 */
  register int j;
  temp2=0.0;
  for(i=0;i<4;i++){
    for(j=0;j<DIMF;j++){
        temp1 = a->d[i].c[j].real * b->d[i].c[j].real; temp2 += temp1;
        temp1 = a->d[i].c[j].imag * b->d[i].c[j].imag; temp2 += temp1;
    }
  }
#endif
  return(temp2);

#else

#ifndef NATIVEDOUBLE
  register double ar,ai,br,bi,ss;
#else
  register Real ar,ai,br,bi,ss;
#endif

  ar=a->d[0].c[0].real;  ai=a->d[0].c[0].imag;
  br=b->d[0].c[0].real;  bi=b->d[0].c[0].imag;
  ss = ar*br + ai*bi;
  ar=a->d[0].c[1].real;  ai=a->d[0].c[1].imag;
  br=b->d[0].c[1].real;  bi=b->d[0].c[1].imag;
  ss += ar*br + ai*bi;
  ar=a->d[0].c[2].real;  ai=a->d[0].c[2].imag;
  br=b->d[0].c[2].real;  bi=b->d[0].c[2].imag;
  ss += ar*br + ai*bi;

  ar=a->d[1].c[0].real;  ai=a->d[1].c[0].imag;
  br=b->d[1].c[0].real;  bi=b->d[1].c[0].imag;
  ss += ar*br + ai*bi;
  ar=a->d[1].c[1].real;  ai=a->d[1].c[1].imag;
  br=b->d[1].c[1].real;  bi=b->d[1].c[1].imag;
  ss += ar*br + ai*bi;
  ar=a->d[1].c[2].real;  ai=a->d[1].c[2].imag;
  br=b->d[1].c[2].real;  bi=b->d[1].c[2].imag;
  ss += ar*br + ai*bi;

  ar=a->d[2].c[0].real;  ai=a->d[2].c[0].imag;
  br=b->d[2].c[0].real;  bi=b->d[2].c[0].imag;
  ss += ar*br + ai*bi;
  ar=a->d[2].c[1].real;  ai=a->d[2].c[1].imag;
  br=b->d[2].c[1].real;  bi=b->d[2].c[1].imag;
  ss += ar*br + ai*bi;
  ar=a->d[2].c[2].real;  ai=a->d[2].c[2].imag;
  br=b->d[2].c[2].real;  bi=b->d[2].c[2].imag;
  ss += ar*br + ai*bi;

  ar=a->d[3].c[0].real;  ai=a->d[3].c[0].imag;
  br=b->d[3].c[0].real;  bi=b->d[3].c[0].imag;
  ss += ar*br + ai*bi;
  ar=a->d[3].c[1].real;  ai=a->d[3].c[1].imag;
  br=b->d[3].c[1].real;  bi=b->d[3].c[1].imag;
  ss += ar*br + ai*bi;
  ar=a->d[3].c[2].real;  ai=a->d[3].c[2].imag;
  br=b->d[3].c[2].real;  bi=b->d[3].c[2].imag;
  ss += ar*br + ai*bi;

  return(ss);

#endif
}
