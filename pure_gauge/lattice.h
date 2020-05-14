// -----------------------------------------------------------------
// Define global scalars and fields in the lattice
#ifndef _LATTICE_H
#define _LATTICE_H

#include "defines.h"
#include "../include/macros.h"  /* For MAXFILENAME, EVENFIRST */
#include "../include/random.h"  /* For double_prn */
#include "../include/io_lat.h"  /* For gauge_file */
#include "../include/su3.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// The lattice is an array of this site struct
typedef struct {
  short x, y, z, t;   // Coordinates of this site
  char parity;        // Is it even or odd?
  int index;          // Index in the array

#ifdef SITERAND
  /* The state information for a random number generator */
  double_prn site_prn;
  /* align to double word boundary (kludge for Intel compiler) */
  int space1;
#endif

    /* Now come the physical fields, program dependent */
  /* gauge field */
  su3_matrix link[4];
  su3_matrix tempmat1,staple;
  /* temporary matrices */
  su3_vector tempvec;  /* for gaugefix */
#endif

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
EXTERN int iseed;    /* random number seed */
EXTERN int warms,trajecs,steps,stepsQ,propinterval;
EXTERN Real beta;
EXTERN Real epsilon;
EXTERN char startfile[MAXFILENAME],savefile[MAXFILENAME];
EXTERN double g_ssplaq, g_stplaq;
EXTERN double_complex linktrsum;
EXTERN u_int32type nersc_checksum;
EXTERN char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable **/
EXTERN int startflag;  /* beginning lattice: CONTINUE, RELOAD, FRESH */
EXTERN int fixflag;  /* gauge fix: COULOMB_GAUGE_FIX, NO_GAUGE_FIX */
EXTERN int saveflag; /* do with lattice: 1=save; */
EXTERN int total_iters;

// Some of these global variables are node dependent
// They are set in "make_lattice()"
EXTERN int sites_on_node;    /* number of sites on this node */
EXTERN int even_sites_on_node; /* number of even sites on this node */
EXTERN int odd_sites_on_node;  /* number of odd sites on this node */
EXTERN int number_of_nodes;  /* number of nodes in use */
EXTERN int this_node;    /* node number of this node */

EXTERN gauge_file *startlat_p;

/* Each node maintains a structure with the pseudorandom number
   generator state */
EXTERN double_prn node_prn ;

// The lattice is a single global variable
// (actually this is the part of the lattice on this node)
EXTERN Real boundary_phase[4];
EXTERN site *lattice;

// Vectors for addressing
// Generic pointers, for gather routines
#define N_POINTERS 8
EXTERN char **gen_pt[N_POINTERS];
#endif
// -----------------------------------------------------------------
