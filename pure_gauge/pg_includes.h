// -----------------------------------------------------------------
// Include files for pure gauge (Wilson action) evolution
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
#include "../include/dirs.h"
#include "../include/field_alloc.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Prototypes for functions in high level code
int setup();
int readin(int prompt);

// Evolution stuff
int update(void);
void update_h(Real eps);
void update_u(Real eps);
void relax(int NumStp);
void monte(int NumStp);
void dsdu_qhb(int dir1, int parity);    // Gauge force for quasi-heat bath
void gauge_field_copy(field_offset src, field_offset dest);
// -----------------------------------------------------------------
