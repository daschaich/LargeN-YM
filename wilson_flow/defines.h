// -----------------------------------------------------------------
// Compiler macros common to all targets in this application
#ifndef _DEFINES_H
#define _DEFINES_H

// No random numbers
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// nHYP stuff
#define NHYP_DEBUG

// Three smearing levels gives nHYP smearing.  This is the maximum
// One smearing level gives something more like stout smearing
// Three alpha_smear parameters must always be supplied
// If SMEAR_LEVEL < 3, extra alpha_smear parameters will be ignored
#define SMEAR_LEVEL 3

// IR regulator:  Q = Omega^dag Omega + IR_STAB
// This slightly changes the definition of the nHYP link. Fine.
// Since we're adding a constant, any derivative of Q is unchanged.
// Used in block_nhyp.c
#define IR_STAB 1e-6

// Neighborhoods of 0 where we use approximate R and S for u0, u1, p
// Bypass R / (-S)^{3 / 2} = 0 / 0 in calculation of Q^{-1 / 2}
// Used in nhyp.c
#if NCOL == 3
  #define EPS_SQ 1e-5
#endif
#if NCOL == 4
  #define EPS_SQ_4 1e-4  // Control fourfold degeneracy
  #define EPS_SQ_3 1e-5  // Control threefold degeneracy
  #ifdef NHYP_JACOBI
    #define MAX_JACOBI_ITERS 1000
    #define TOL_JACOBI 1e-8
  #endif
#endif

#ifdef NHYP_DEBUG
  #if NCOL == 4
    #define TOL_NHYP 1e-5          // Tolerance for warnings
  #endif
#endif

#endif
// -----------------------------------------------------------------
