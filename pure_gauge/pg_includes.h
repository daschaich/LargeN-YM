// -----------------------------------------------------------------
// Include files for pure gauge (Wilson action) evolution
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>             // For strlen
#include <math.h>
#include <time.h>
#include <math.h>
#include "../include/config.h"  // Keep this first
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/macros.h"
#include "lattice.h"
#include "../include/comdefs.h"
#include "../include/io_lat.h"
#include "../include/generic.h"
#include "../include/dirs.h"
#include "../include/field_alloc.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Prototypes for functions in high level code
int setup();
int readin(int prompt);

// Generic update routine switches between ora and HMC
void update();

#ifndef HMC
// Over-relaxed quasi-heatbath stuff
void relax();
void dsdu_qhb(int dir, int parity);    // Gauge force for quasi-heatbath
void monte();
#else
// HMC stuff
double U_action();
double action();
double action_HMC();
void update_u(Real eps);
double update_h(Real eps);
int update_hmc();
#endif

// LLR stuff
#ifdef LLR
double energy(double *ss_plaq, double *st_plaq);
void findEint();

void updateconst_e();
#ifdef HMC
  // HMC updates with gaussian window
double update_h_const();
int update_hmc_const();
#else
  // Over-relaxed quasi-heatbath updates with hard constraints
void monteconst_e();
#endif
#endif
// -----------------------------------------------------------------
