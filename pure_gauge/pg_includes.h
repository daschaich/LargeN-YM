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

// Over-relaxed quasi-heatbath stuff
void relax();
void dsdu_qhb(int dir1, int parity);    // Gauge force for quasi-heat bath
void monte();
void update();

// HMC stuff
 void gauge_field_copy(field_offset src, field_offset dest);
double action();
int update_hmc();
void update_u();
double update_h();

#ifdef LLR
// Impose constraint that energy remains within interval under consideration
void updateconst_e(double Eint, double a);
void monteconst_e(double Eint, double a);
#endif
// -----------------------------------------------------------------
