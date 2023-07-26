// -----------------------------------------------------------------
// Define global scalars and fields in the lattice
#ifndef _LATTICE_H
#define _LATTICE_H

#include "defines.h"
#include "../include/macros.h"    // For MAXFILENAME
#include "../include/io_lat.h"    // For gauge_file
#include "../include/random.h"    // For double_prn
#include "../include/su3.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// The lattice is an array of this site struct
typedef struct {
  short x, y, z, t;   // Coordinates of this site
  char parity;        // Is it even or odd?
  int index;          // Index in the array

#ifdef SITERAND
  // The state information for a random number generator
  double_prn site_prn;
  // Align to double word boundary (kludge for Intel compiler)
  int space1;
#endif

  // Gauge links
  matrix link[4];
#ifdef HMC
  matrix old_link[4];  // For accept/reject

  // Antihermitian momentum matrices in each direction
  anti_hermitmat mom[4];
#endif

  // Temporary storage
  matrix tempmat, staple;    // TODO: Convert these to fields
} site;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Definition of global variables
#ifdef CONTROL
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN int nx, ny, nz, nt;  // Lattice dimensions
EXTERN int volume;          // Volume of lattice
EXTERN int iseed;           // Random number seed
EXTERN int warms, trajecs;  // Common stuff
EXTERN int measinterval;

EXTERN int ora_steps, qhb_steps;    // ORA stuff
#ifndef HMC
//EXTERN int ora_steps, qhb_steps;    // ORA stuff
#else
EXTERN int hmc_steps;  // HMC stuff
EXTERN Real traj_length;  // HMC stuff
EXTERN Real fnorm, max_f;           // Force monitoring
#endif

EXTERN Real beta;
EXTERN Real one_ov_N, one_ov_vol;
EXTERN char startfile[MAXFILENAME], savefile[MAXFILENAME];
EXTERN double g_ssplaq, g_stplaq;
EXTERN double_complex linktr;
EXTERN u_int32type nersc_checksum;
EXTERN char stringLFN[MAXFILENAME];   // ILDG LFN if applicable
EXTERN int startflag; // Beginning lattice: CONTINUE, RELOAD, FRESH
EXTERN int fixflag;   // Either NO_GAUGE_FIX or COULOMB_GAUGE_FIX
EXTERN int saveflag;  // 1 if we will save the lattice
EXTERN int total_iters;

#ifdef LLR
// LLR parameters
EXTERN int ait, accept, reject, constrained, Njacknife;
EXTERN Real a, Emin, Emax, delta, deltaSq;
#endif

// Some of these global variables are node dependent
// They are set in "make_lattice()"
EXTERN int sites_on_node;       // Number of sites on this node
EXTERN int even_sites_on_node;  // Number of even sites on this node
EXTERN int odd_sites_on_node;   // Number of odd sites on this node
EXTERN int number_of_nodes;     // Number of nodes in use
EXTERN int this_node;           // Node number of this node

// Each node maintains a structure with the pseudorandom number
// generator state
EXTERN double_prn node_prn;

// Temporary fields
EXTERN matrix *tempmat, *tempmat2;

// Some more arrays to be used by LAPACK
// in reunitarization (in generic directory)
EXTERN double *Rwork, *eigs, *store, *work, *junk, *left, *right;

// The lattice is a single global variable
// (actually this is the part of the lattice on this node)
EXTERN Real boundary_phase[4];
EXTERN site *lattice;

EXTERN gauge_file *startlat_p;

// Vectors for addressing
// Generic pointers, for gather routines
#define N_POINTERS 8
EXTERN char **gen_pt[N_POINTERS];
#endif
// -----------------------------------------------------------------
