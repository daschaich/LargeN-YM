/*******************  trace_su3.c  (in su3.a) ***************************
*									*
* complex trace_su3(a) su3_matrix *a;					*
* return complex trace of an SU3 matrix 				*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

/* Complex trace of an SU3 matrix */
complex trace_su3( su3_matrix *a ) {
register complex t1,t2;
#if (DIMF<=6)
    CADD(a->e[0][0],a->e[1][1],t1);
#if (DIMF==2)
    return(t1);
#else
    CADD(t1,a->e[2][2],t2);
#if (DIMF==3)
    return(t2);
#else
    CADD(t2,a->e[3][3],t1);
#if (DIMF==4)
    return(t1);
#else
    CADD(t1,a->e[4][4],t2);
#if (DIMF==5)
    return(t2);
#else
    CADD(t2,a->e[5][5],t1);
    return(t1);
#endif
#endif
#endif
#endif
#else  /* DIMF>6 */
    int i;
    CADD(a->e[0][0],a->e[1][1],t1);
    for (i=2;i<DIMF;i++){
        CADD(t1,a->e[i][i],t1);
    }
    return(t1);
#endif
}
