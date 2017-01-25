/*********************  subvec.c  (in su3.a) ****************************
*									*
* void sub_vector(a,b,c) vector *a,*b,*c;			*
* subtract su3 vectors:  C <-  A - B 					*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

/* subtract su3 vectors */
void sub_vector( vector *a, vector *b, vector *c ){
register int i;
    for(i=0;i<DIMF;i++){
	CSUB( a->c[i], b->c[i], c->c[i] );
    }
}
