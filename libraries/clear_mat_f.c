/********************  clear_mat_f.c  (in su3.a) ********************
*
*void clear_su3mat_f( su3_matrix_f *dest )
*  clear an SU3 matrix
* dest  <-  zero_matrix
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void clear_su3mat_f( su3_matrix_f *dest ){
register int i,j;
    for(i=0;i<NCOL;i++)for(j=0;j<NCOL;j++){
	dest->e[i][j].real = dest->e[i][j].imag = 0.0;
    }
}
