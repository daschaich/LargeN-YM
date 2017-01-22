// -----------------------------------------------------------------
// Include files for dynamical Wilson-clover HMC
#include <stdio.h>
#include <stdlib.h>
#include <string.h>             // For strlen
#include <math.h>
#include "../include/config.h"  // Keep this first
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/macros.h"
#include "lattice.h"
#include "../include/comdefs.h"
#include "../include/io_lat.h"
#include "../include/generic.h"
#include "../include/generic_wilson.h"
#include "../include/generic_clover.h"
#include "../include/generic_nhyp.h"
#include "../include/dirs.h"
#include "../include/field_alloc.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Prototypes for functions in high-level code
int setup();
int readin(int prompt);

// CG stuff
void fermion_op(wilson_vector *src, wilson_vector *dest,
                int sign, int parity);
int congrad(int level, Real mshift, int parity);

// Measurement stuff
double action();
int f_meas();
#ifdef SPECTRUM
int spectrum();
#endif

// Evolution stuff
void grsource();
void checkmul(int level, Real mshift);
int update();
void update_u(Real eps);
void update_h(Real eps);
void udadu_mat_mu_nu(int mu, int nu);
void udadu_mu_nu(wilson_vector *lsrc, wilson_vector *rsrc, int mu, int nu);
void mult_sigma_mu_nu(wilson_vector *src, wilson_vector *dest,
                      int mu, int nu);
void chain_rule(su3_matrix_f *sigmaf, su3_matrix *sigma,
                su3_matrix_f *gaugelinkf);
void apply_bc(su3_matrix_f *sigmaf, int dir, int t);
double gauge_force(Real eps);
double fermion_force(Real eps1, Real eps2);

// nHYP force stuff
void nhyp_force1();
// -----------------------------------------------------------------
