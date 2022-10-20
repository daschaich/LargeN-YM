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
  matrix_f linkf[4];
#ifdef HMC
  matrix_f old_linkf[4];  // For accept/reject

  // Antihermitian momentum matrices in each direction
  anti_hermitmat mom[4];
#endif

  // Temporary storage
  matrix_f tempmat, staple;    // TODO: Replace with tempmatf and tempmatf2
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
#ifndef HMC
EXTERN int ora_steps, qhb_steps;    // ORA stuff
#else
EXTERN int hmc_steps, traj_length;  // HMC stuff
EXTERN Real fnorm, max_f;           // Force monitoring
#endif

EXTERN Real beta, a;        // Fix a = 1 without LLR
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
EXTERN int ait, accept, reject;
EXTERN double Emin, Emax, delta, deltaSq;
#endif

// Some of these global variables are node dependent
// They are set in "make_lattice()"
EXTERN int sites_on_node;       // Number of sites on this node
EXTERN int even_sites_on_node;  // Number of even sites on this node
EXTERN int odd_sites_on_node;   // Number of odd sites on this node
EXTERN int number_of_nodes;     // Number of nodes in use
EXTERN int this_node;           // Node number of this node

EXTERN gauge_file *startlat_p;

// Each node maintains a structure with the pseudorandom number
// generator state
EXTERN double_prn node_prn;

// The lattice is a single global variable
// (actually this is the part of the lattice on this node)
EXTERN Real boundary_phase[4];
EXTERN site *lattice;

// Vectors for addressing
// Generic pointers, for gather routines
#define N_POINTERS 8
EXTERN char **gen_pt[N_POINTERS];
#endif

// Temporary fields for field_strength
EXTERN matrix_f *tempmatf, *tempmatf2;
// -----------------------------------------------------------------
