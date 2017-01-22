#ifndef _GENERIC_QUARK_TYPES_H
#define _GENERIC_QUARK_TYPES_H

#include "../include/macros.h"
#include "../include/precision.h"

// Struct to hold quark inversion parameters
typedef struct {
  int min;            // Minimum number of iterations
  int max;            // Maximum number of iterations per restart
  int nrestart;       // Maximum restarts
  int parity;         // EVEN, ODD or EVENANDODD
  int start_flag;     // Zero out initial guess if 0, otherwise use psi[0]
  int nsrc;           // Number of source vectors
  Real resid;         // Target residual, normalized as
                      //   sqrt(|r|^2) / sqrt(|src_e|^2)
  Real size_r;        // Final residual
  int converged;      // 1 of converged, otherwise 0
} quark_invert_control;

/* Struct defining parameters of Dirac matrix for clover inversion */
/* To be passed through to inverter */
typedef struct {
  Real Kappa;        // Hopping parameter
  Real Clov_c;       // Clover coefficient
  Real U0;           // Tadpole coefficient
} dirac_clover_param;

#endif
