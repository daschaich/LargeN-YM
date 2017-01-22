// -----------------------------------------------------------------
// Macros and declarations for generic clover routines
#ifndef _GENERIC_CLOVER_H
#define _GENERIC_CLOVER_H

#include "../include/generic_quark_types.h"
#include "../include/macros.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void f_mu_nu(su3_matrix f_mn[],int mu,int nu);

typedef struct { complex tr[2][DIMF*(2*DIMF-1)]; } triangular;
typedef struct { Real di[2][2*DIMF]; } diagonal;
typedef struct {
  triangular *clov;
  diagonal *clov_diag;
} clover;

/* make_clov.c routines for any clover term */
clover *create_clov(void);
void compute_clov(clover *my_clov, Real Clov_c);
double compute_clovinv(clover *my_clov, int parity);

// Recast src and dest as wilson_block_vector
void mult_this_ldu(clover *my_clov, wilson_vector *src,
                   wilson_vector *dest, int parity);
void free_this_clov(clover *my_clov);

/* make_clov.c routines for single clover term */
void make_clov(Real Clov_c);
double make_clovinv(int parity);
void tr_sigma_ldu_mu_nu(int mu, int nu);
void free_clov();
void mult_ldu(wilson_vector *src, wilson_vector *dest, int parity);
#endif
// -----------------------------------------------------------------
