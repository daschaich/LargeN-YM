/******************  realtr.c  (in su3.a) *******************************
*									*
* Real realtrace_su3( su3_matrix *a,*b)				*
* return Re( Tr( A_adjoint*B )  					*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

Real realtrace_su3(  su3_matrix *a, su3_matrix *b ){
register int i,j;
register Real sum;
    for(sum=0.0,i=0;i<DIMF;i++)for(j=0;j<DIMF;j++)
	sum+= a->e[i][j].real*b->e[i][j].real + a->e[i][j].imag*b->e[i][j].imag;
    return(sum);
}
