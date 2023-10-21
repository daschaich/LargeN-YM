// -----------------------------------------------------------------
// Compiler macros common to all targets in this application
#ifndef _DEFINES_H
#define _DEFINES_H

#define PURE_GAUGE            // No fermions
#define SITERAND              // Use site-based random number generators
#define GAUGE_FIX_TOL 1.0e-7  // For gauge fixing

#ifdef LLR
// LLR stuff
#define NRiters 30
#define a_cut 200.0
#define FIND_MAX 2000
#endif

#endif
// -----------------------------------------------------------------
