/*******************  submat_f.c  (in su3.a) ******************************
*									*
* void sub_su3_matrix_f(a,b,c) su3_matrix_f *a,*b,*c;			*
* subtract su3 matrices:  C  <- A - B 					*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

/* subtract su3 matrices */
void sub_su3_matrix_f( su3_matrix_f *a, su3_matrix_f *b, su3_matrix_f *c ) {
register int i,j;
    for(i=0;i<NCOL;i++)for(j=0;j<NCOL;j++){
	CSUB( a->e[i][j], b->e[i][j], c->e[i][j] );
    }
}
