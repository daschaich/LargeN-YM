/*******************  clearvec_f.c  (in su3.a) *****************************
*									*
*  void clearvec_f( su3_vector_f *vec )					*
*  clear a NCOL element complex vector					*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void clearvec_f( su3_vector_f *v ){
     v->c[0].real = v->c[0].imag = 0.0;
     v->c[1].real = v->c[1].imag = 0.0;
#if (NCOL>2)
     v->c[2].real = v->c[2].imag = 0.0;
#if (NCOL>3)
     v->c[3].real = v->c[3].imag = 0.0;
#endif
#endif
}
