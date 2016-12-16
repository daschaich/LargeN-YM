/*******************  trace_su3_f.c  (in su3.a) ***************************
*									*
* complex trace_su3_f(a) su3_matrix_f *a;					*
* return complex trace of an SU3 matrix 				*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

/* Complex trace of an SU3 matrix */
complex trace_su3_f( su3_matrix_f *a ) {
register complex t1,t2;
    CADD(a->e[0][0],a->e[1][1],t1);
#if (NCOL==2)
    return(t1);
#else
    CADD(t1,a->e[2][2],t2);
#if (NCOL==3)
    return(t2);
#else
    CADD(t2,a->e[3][3],t1);
    return(t1);
#endif
#endif
}
