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

#if defined(ORA) || defined(LLR)
// Over-relaxed quasi-heatbath stuff
// Used to find energy interval when doing LLR with HMC
void dsdu_qhb(int dir, int parity);    // Gauge force for quasi-heatbath
void relax();
void monte();
void update_ora();
#endif

// HMC stuff, including energy interval info for LLR
#ifdef HMC
double action(double E_min);
void update_u(Real eps);
double update_h(Real eps, double E_min);
void hmc_traj(double E_min);
void update_hmc(double E_min);
#endif

// LLR stuff
#ifdef LLR
double gauge_action();
void findEint(double E_min);
#endif
// -----------------------------------------------------------------
