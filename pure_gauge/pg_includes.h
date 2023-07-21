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

void relax();
void monte();
void monteconst_e();
void frankensteinconst_e(double E_min,int node_counter);
void dsdu_qhb(int dir, int parity);    // Gauge force for quasi-heatbath
#ifndef HMC
// Over-relaxed quasi-heatbath stuff
#else
// HMC stuff
double action(double E_min);
void update_u(Real eps);
double update_h(Real eps,double E_min);
void update_hmc(double E_min);
#endif

// Need gauge_action() for LLR with or without HMC
#if defined(HMC) || defined(LLR)
double gauge_action();
#endif

// LLR stuff
#ifdef LLR
void findEint(double E_min);

void updateconst_e(double E_min);
#ifndef HMC
  // Over-relaxed quasi-heatbath updates with hard constraints

#endif
#endif
// -----------------------------------------------------------------
