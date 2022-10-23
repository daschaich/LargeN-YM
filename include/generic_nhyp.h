// -----------------------------------------------------------------
// Includes and definitions for nHYP smearing
#ifndef _GENERIC_NHYP_H
#define _GENERIC_NHYP_H
#include "../include/su3.h"
#include "../include/comdefs.h"
#include "../include/macros.h"
#include "../include/generic_quark_types.h"

void block_nhyp();

// Code related to b[][] depends on NCOL
#ifndef NHYP_DEBUG
void compute_fhb(matrix *Q, Real *f, Real b[NCOL][NCOL], int compute_b);
#else
void compute_fhb(matrix *Omega, matrix *Q, Real *f,
                 Real b[NCOL][NCOL], int compute_b);
#endif

#endif // _GENERIC_NHYP_H
// -----------------------------------------------------------------
