// -----------------------------------------------------------------
// Define global scalars and fields in the lattice
#ifndef _LATTICE_H
#define _LATTICE_H

#include "defines.h"
#include "../include/generic_quark_types.h"
#include "../include/macros.h"    // For MAXFILENAME
#include "../include/io_lat.h"    // For gauge_file
#include "../include/random.h"    // For double_prn
#include "../include/su3.h"

#ifdef NHYP_JACOBI
#include "../include/jacobi.h"
#endif
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

  // Gauge field
  matrix link[4];
  matrix_f linkf[4];

  // Program-dependent fields
#ifdef HMC_ALGORITHM
  matrix_f old_linkf[4];    // For accept/reject
#endif

  // Antihermitian momentum matrices in each direction
  anti_hermitmat mom[4];

#ifdef SPECTRUM
  wilson_matrix quark_propagator;
  wilson_matrix rotated_propagator;   // For clover-rotated operators
#endif

  // A bit wasteful of space
  matrix Force[4];
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
EXTERN int warms, trajecs, niter, nrestart, propinterval, nflavors;
EXTERN float traj_length;

// Global Hasenbusch variables
EXTERN int num_masses;    // Maximum number of masses, <= MAX_MASSES
EXTERN Real shift;
EXTERN complex ishift, mishift;
EXTERN int nsteps[MAX_MASSES + 1];
EXTERN Real gnorm, fnorm[2], max_gf, max_ff[2];

// Global action and evolution variables
EXTERN Real beta, beta_frep, kappa, mkappa, mkappaSq, clov_c, u0, CKU0;
EXTERN Real rsqmin, rsqprop, one_ov_N, one_ov_vol;
EXTERN double returntrlogA, g_ssplaq, g_stplaq;
EXTERN double_complex linktr;
EXTERN u_int32type nersc_checksum;
EXTERN char stringLFN[MAXFILENAME];   // ILDG LFN if applicable
EXTERN char startfile[MAXFILENAME], savefile[MAXFILENAME];
EXTERN int startflag; // Beginning lattice: CONTINUE, RELOAD, FRESH
EXTERN int fixflag;   // Either NO_GAUGE_FIX or COULOMB_GAUGE_FIX
EXTERN int saveflag;  // 1 if we will save the lattice
EXTERN int total_iters;
/* boundary-flip switch, made consistent for higher-rep code */
EXTERN int current_boundary;

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

EXTERN gauge_file *startlat_p;
EXTERN gauge_file *savelat_p;

// The lattice is a single global variable
// (actually this is the part of the lattice on this node)
EXTERN site *lattice;

// Vectors for addressing
// Generic pointers, for gather routines
// We need 8 mostly (since force_nhyp is streamlined)
#define N_POINTERS 8
EXTERN char **gen_pt[N_POINTERS];

EXTERN matrix_f *gauge_field[4];
EXTERN matrix_f *gauge_field_thin[4];

// CG stuff
EXTERN wilson_vector *g_rand;           // Gaussian random vector
EXTERN wilson_vector *chi[MAX_MASSES];  // Gaussian random source vector
EXTERN wilson_vector *psi[MAX_MASSES];  // Solution vector = Minverse * chi
EXTERN wilson_vector *p;                // CG change vector
EXTERN wilson_vector *mp;               // CG temporary vector
EXTERN wilson_vector *r;                // CG residual vector
EXTERN wilson_vector *tempwvec;         // Another temporary CG vector
#ifdef PHI_ALGORITHM
EXTERN wilson_vector *old_psi[MAX_MASSES];  // For predicting next psi
#endif
EXTERN half_wilson_vector *htmp[8];     // Temporaries for dslash_w_field

// Clover stuff
EXTERN matrix *f_mn;

// nHYP stuff
EXTERN Real alpha_smear[3];
EXTERN matrix_f *hyplink1[4][4];
EXTERN matrix_f *hyplink2[4][4];
EXTERN matrix_f *Sigma[4];
EXTERN matrix_f *SigmaH[4];
EXTERN matrix_f *SigmaH2[4][4];

EXTERN matrix_f *Staple1[4][4];
EXTERN matrix_f *Staple2[4][4];
EXTERN matrix_f *Staple3[4];

EXTERN matrix_f *LambdaU[4];
EXTERN matrix_f *Lambda1[4];
EXTERN matrix_f *Lambda2[4];

// Temporary fundamental and irrep matrices
EXTERN matrix_f *tempmatf, *tempmatf2, *staplef;
EXTERN matrix *tempmat, *tempmat2, *staple;

#ifdef NHYP_JACOBI
EXTERN Matrix Qj, Vj;
#define JACOBI_HIST_MAX 10
EXTERN int jacobi_hist[JACOBI_HIST_MAX];
#endif

// Up to 20 concurrent timers for timing
#ifdef TIMING
EXTERN double tmptime[20];
EXTERN double time_dcongrad;      // From congrad.  Use tmptime[0]
EXTERN double time_fermion_force; // From update.  Use tmptime[1]
EXTERN double time_fermion_rep;   // From update.  Use tmptime[2]
EXTERN double time_block_nhyp;    // From block_nhyp.  Use tmptime[3]
EXTERN double time_compute_fhb;   // From compute_fhb.  Use tmptime[4]
EXTERN double time_jacobi;        // Not currently using tmptime[6]
#endif

#endif
// -----------------------------------------------------------------
