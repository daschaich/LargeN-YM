/******************  su3_adjoint_f.c  (in su3.a) **************************
*									*
* void su3_adjoint_f( su3_matrix_f *a, su3_matrix_f *b )		*
* B  <- A_adjoint,  adjoint of an SU3 matrix 				*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

/* adjoint of an SU3 matrix */
void su3_adjoint_f( su3_matrix_f *a, su3_matrix_f *b ){
register int i,j;
    for(i=0;i<NCOL;i++)for(j=0;j<NCOL;j++){
	CONJG( a->e[j][i], b->e[i][j] );
    }
}
