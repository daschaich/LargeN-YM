/********************  copy_wvec.c  (in su3.a) ********************
*
*void copy_wvec( wilson_vector *src,*dest )
*  copy a Wilson vector
* dest  <-  src
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

/* IS THIS CORRECT ?? */

void copy_wvec( wilson_vector *src, wilson_vector *dest ){
    *dest = *src;	/* hardly worth a function */
}
