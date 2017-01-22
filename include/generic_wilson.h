// -----------------------------------------------------------------
// Macros and declarations for generic Wilson routines
#ifndef _GENERIC_WILSON_H
#define _GENERIC_WILSON_H

#include "../include/su3.h"
#include "../include/comdefs.h"
#include "../include/macros.h"
#include "../include/generic_quark_types.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void boundary_flip(int sign);

void dslash(wilson_vector *src, wilson_vector *dest, int isign, int parity);
void dslash_special(wilson_vector *src, wilson_vector *dest,
                    int isign, int parity, msg_tag **tag, int is_started);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Higher-rep routines
void fermion_rep();
void make_fermion_rep_matrix(su3_matrix_f *a, su3_matrix *b);
#endif
// -----------------------------------------------------------------
