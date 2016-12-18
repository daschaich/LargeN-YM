// -----------------------------------------------------------------
// Include files for dynamical Wilson-clover HMC
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../include/config.h"  /* Keep this first */
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

void gauge_field_copy(field_offset src,field_offset dest);
void gauge_field_copy_f(field_offset src,field_offset dest);

// CG stuff
int congrad_cl_m(int niter, Real rsqmin, Real *final_rsq_ptr,
                 field_offset src, field_offset dest, Real mshift);

void gauge_action(double *result);
double d_action();
int f_measure_cl();
int pcac_t();
int pcac_x();
int w_spectrum_cl();
int t_props_cl();
int s_props_cl();
void make_loop_table2();
void path(int *dir,int *sign,int length);
void single_action(int dir, Real *coeff);
void udadu_mu_nu(field_offset lsrc, field_offset rsrc, field_offset mat,
     int mu, int nu, int parity);
void udadu_mat_mu_nu(field_offset matsrc, field_offset matdest,
         int mu, int nu);
void checkmul(field_offset chi, field_offset psi, Real mshift);

void prepare_vecs(int level);
void chain_rule(su3_matrix_f *sigmaf, su3_matrix *sigma,
                su3_matrix_f *gaugelinkf);
void apply_bc(su3_matrix_f *sigmaf, int dir, int t);
void mult_sigma_mu_nu(wilson_vector *src, wilson_vector *dest,
          int mu, int nu);

// Evolution stuff
void grsource_w();
int update();
void update_h(Real eps);
void update_u(Real eps);
double gauge_force(Real eps);
double fermion_force(Real eps1, Real eps2);

// nHYP force stuff
void stout_force1();
// -----------------------------------------------------------------
