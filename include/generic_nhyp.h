#ifndef _GENERIC_NHYP_H
#define _GENERIC_NHYP_H
/************************ generic_nhyp.h ******************************
*									*
*  Macros and declarations for generic_nhyp routines                    *
*  This header is for codes that call generic_nhyp routines             *
*/

#include "../include/su3.h"
#include "../include/comdefs.h"
#include "../include/macros.h"
#include "../include/generic_quark_types.h"

/* code related to b[][] is specific to SU(2), SU(3)  */

#ifndef NHYP_DEBUG
void compute_fhb( su3_matrix_f *Q, Real *f, Real b[NCOL][NCOL], int compute_b );
#else
void compute_fhb( su3_matrix_f *Omega, su3_matrix_f *Q,
                  Real *f, Real b[NCOL][NCOL], int compute_b );
#endif
void block_nhyp();


#endif /* _GENERIC_NHYP_H */
