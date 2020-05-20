// -----------------------------------------------------------------
// Include files for SU(N) Wilson flow measurements
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
#include "../include/generic_nhyp.h"
#include "../include/dirs.h"
#include "../include/field_alloc.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Prototypes for functions in high level code
int setup();
int readin(int prompt);
void wflow();
void stout_step_rk();
void staple(matrix_f *stp[NDIMS]);

// Measurement stuff
void meas(Real t);
void make_loop_table();
void gauge_loops(double *result);
void path(int *dir, int *sign, int length);
// -----------------------------------------------------------------
