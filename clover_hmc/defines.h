// -----------------------------------------------------------------
// Compiler macros common to all targets in this application
#ifndef _DEFINES_H
#define _DEFINES_H

#define SITERAND   /* Use site-based random number generators */
#define GAUGE_FIX_TOL 1.0e-7 /* For gauge fixing */

#define TIMING
//#define CG_DEBUG

/* number of CG restarts during update                              */
#define CONGRAD_RESTART 10


// Maximal number of masses used in the Hasenbusch preconditioning
#define MAX_MASSES 2

/* flags for pcac_x                                                 */
/* measure everything also with PBC+APBC                            */
#define D_TOO
/* measure <V1(0)V1(x)>                                             */
#define V_TOO

/* in pcac_sf, measure also correlators with UP boundary */
/* #define UP_TOO */

/**************************** SF stuff ******************************/

#ifdef SF
/*  x-direction background field: 1 turns on, 0 turns off           */
#define BACKGROUND_FIELD_X 1
#endif

/*********************** fat links stuff ****************************/

#ifdef NHYP

/*** number of smearing levels ***/

/* Set the number of smearing levels. Legal values are 1, 2, or 3
   (true NHYP is 3, "stout" is 1).
   NOTE: to save coding, the three alpha_smear parameters must
   always be supplied, either by hard coding below or by reading them from
   the infile.  But if SMEAR_LEVEL=n with n<3, then only the first n values
   of the alpha_smear parameters are used.
*/
#define SMEAR_LEVEL 3

/*** hard coding smearing parameters ***/

/* Comment out the following line
   if you want to read the smearing parameters from the infile.
   Otherwise, they are hard-coded below
*/
/*
   #define HARD_CODE_SMEAR
*/
#ifdef HARD_CODE_SMEAR
#ifdef CONTROL
Real alpha_smear[3]={0.5, 0.5, 0.4};
#else
extern Real alpha_smear[3];
#endif /* CONTROL */
#endif /* HARD_CODE_SMEAR */

/*** constants for numerical stability ***/

/* IR regulator:  Q = Omega^dag Omega + IR_STAB
   This slightly changes the definition of the nhyp link. Fine.
   Since we're adding a constant, any derivative of Q is unchanged.
   routines block_nhyp1,2,3() in block_nhyp.c
            Sigma_update1() in force_nhyp.c
*/
#define IR_STAB 1.0e-6

/* calculation of Q^{-1/2}:  bypass R/(-S)^{3/2}=0/0
   routine compute_fhb() in generic_nhyp/nhyp.c
   neighborhood of 0 where we use approximate R and S for u0, u1, p
*/

#if(NCOL==3)
  #define EPS_SQ 1.0e-5
#endif
#if(NCOL==4)
  #define EPS_SQ_4 1.0e-4  /* control fourfold degeneracy  */
  #define EPS_SQ_3 1.0e-5  /* control threefold degeneracy */
  #ifdef NHYP_JACOBI
    #define MAX_JACOBI_ITERS 1000
    #define TOL_JACOBI 1.e-8
  #endif
#endif

#define NHYP_DEBUG

#ifdef NHYP_DEBUG
#if(NCOL==4)
#define TOL_NHYP 1.0e-5          /* tolerance for warnings     */
#define PRINT_JACOBI_ITERS 20    /* monitor Jacobi iterations  */
/* the following are not used with jacobi */
#define TOL_ACOS 1.0e-5    /* tolerance for acos_input warning */
#define TOL_GMP  1.0e-6    /* monitor jacobi vs gmp            */
#endif
#endif

#endif /* NHYP*/
/******************** END fat links stuff ***************************/
/* This is the lambda for the Omelyan integrator */
#define INT_LAMBDA 0.193
/* This is 2*lambda */
#define INT_LAMBDA_CONT 0.386
/* This is 1-2*lambda */
#define INT_LAMBDA_MID 0.614


#endif /* _DEFINES_H */
