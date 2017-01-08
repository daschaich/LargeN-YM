#ifndef _LATTICE_H
#define _LATTICE_H
/****************************** lattice.h ********************************/

/* defines global scalars and the fields in the lattice. */

#include "defines.h"
#include "../include/generic_quark_types.h"
#include "../include/macros.h"    /* For MAXFILENAME */
#include "../include/io_lat.h"    /* For gauge_file */
#include "../include/random.h"    /* For double_prn */
#include "../include/su3.h"

#ifdef NHYP_JACOBI
#include "../include/jacobi.h"
#endif

// The lattice is an array of sites
typedef struct {
    /* The first part is standard to all programs */
  /* coordinates of this site */
  short x,y,z,t;
  /* is it even or odd? */
  char parity;
  /* my index in the array */
  int index;
#ifdef SITERAND
  /* The state information for a random number generator */
  double_prn site_prn;
  /* align to double word boundary (kludge for Intel compiler) */
  int space1;
#endif

    /* Now come the physical fields, program dependent */
  /* gauge field */
  su3_matrix link[4];
  su3_matrix_f linkf[4];
#ifdef HMC_ALGORITHM
  su3_matrix_f old_linkf[4];    // For accept/reject
#endif

  /* antihermitian momentum matrices in each direction */
  anti_hermitmat mom[4];

  /* Wilson complex vectors */
  wilson_vector g_rand; /* gaussian random vector */
  wilson_vector psi[MAX_MASSES];  /* solution vector */
  wilson_vector chi[MAX_MASSES];  /* source vector */
  wilson_vector p;        /* conjugate gradient change vector */
  wilson_vector mp;       /* another CG vector */
  wilson_vector tmp;      /* another temporary CG vector */
  wilson_vector r;  /* residue */
  wilson_vector invp;     /* used in paca_[tx] */
  wilson_vector inva;     /* used in paca_[tx] */

  /* wilson half vector (temporary used in dslash_w_site) */
  half_wilson_vector htmp[8];
/**half_wilson_vector htmp2[8];**/ /* TEMP FOR TESTING */
#ifdef PHI_ALGORITHM
  wilson_vector old_psi[MAX_MASSES];  // For predicting next psi
#endif
#ifdef SPECTRUM
  wilson_matrix quark_propagator;
      /* For four source spins, three source colors */
  wilson_matrix rotated_propagator;
      /* For clover-rotated operators */
#endif
  /* temporary vectors and matrices */
  su3_matrix tempmat1,tempmat2,staple; /* sometimes these will be
            cast to su3_matrix_f as
            needed for scratch space */

  // A bit wasteful of space
  su3_matrix Force[4];
} site;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Definition of global variables
EXTERN int nx, ny, nz, nt;  // Lattice dimensions
EXTERN int volume;          // Volume of lattice
EXTERN int iseed;           // Random number seed
EXTERN int warms, trajecs, niter, nrestart, propinterval, nflavors;
EXTERN float traj_length;

// Global Hasenbusch variables
EXTERN int num_masses;      // Maximum number of masses, <= MAX_MASSES
EXTERN Real shift;
EXTERN int nsteps[MAX_MASSES + 1];
EXTERN Real max_gf, max_ff[2];

// Global action and evolution variables
EXTERN Real rsqmin, rsqprop, beta, beta_frep, kappa, clov_c, u0;
EXTERN Real epsilon;
EXTERN char startfile[MAXFILENAME],savefile[MAXFILENAME];
EXTERN double g_ssplaq, g_stplaq;
EXTERN double_complex linktr;
EXTERN u_int32type nersc_checksum;
EXTERN char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable **/
EXTERN int startflag;  /* beginning lattice: CONTINUE, RELOAD, FRESH */
EXTERN int fixflag;  /* gauge fix: COULOMB_GAUGE_FIX, NO_GAUGE_FIX */
EXTERN int saveflag; /* do with lattice: 1=save; */
EXTERN int total_iters;
/* boundary-flip switch, made consistent for higher-rep code */
EXTERN int current_boundary;
EXTERN int current_boundary_x;

/* Some of these global variables are node dependent */
/* They are set in "make_lattice()"                  */
EXTERN int sites_on_node;    /* number of sites on this node */
EXTERN int even_sites_on_node; /* number of even sites on this node */
EXTERN int odd_sites_on_node;  /* number of odd sites on this node */
EXTERN int number_of_nodes;  /* number of nodes in use */
EXTERN int this_node;    /* node number of this node */

EXTERN int debugflag;

/* Each node maintains a structure with the pseudorandom number
   generator state */
EXTERN double_prn node_prn ;

/* improved action stuff */
#ifdef IMP

#define nloop 1
#define nreps 1
#ifdef GB
#define nloop_m 4
#define nreps_m 3
#endif
#define max_num 400

/* global defns for general action (check this works...) */
EXTERN int loop_ind[nloop][10], loop_length[nloop];
EXTERN int loop_table[nloop][max_num][10],loop_num[nloop],loop_char[max_num];
EXTERN Real loop_coeff[nloop][nreps];
EXTERN int loop_ch[nloop][max_num],ch;

EXTERN Real loop_term[48][nreps];

#ifdef GB
EXTERN int loop_ind_m[nloop_m][10], loop_length_m[nloop_m];
EXTERN int loop_table_m[nloop_m][max_num][10],loop_num_m[nloop_m],
loop_char_m[max_num];
EXTERN int loop_ch_m[nloop_m][max_num];
#endif

#endif /* IMP */

EXTERN quark_invert_control qic;
EXTERN dirac_clover_param dcp;

EXTERN gauge_file *startlat_p;
EXTERN gauge_file *savelat_p;

/* The lattice is a single global variable - (actually this is the
   part of the lattice on this node) */
EXTERN site *lattice;

/* Vectors for addressing */
/* Generic pointers, for gather routines */
#define N_POINTERS 8  /* Number of generic pointers */
/* NEED 8 WHEN GAUGEFIXING; 10 in original nhyp code */
EXTERN char ** gen_pt[N_POINTERS];

EXTERN su3_matrix_f *gauge_field[4];
EXTERN su3_matrix_f *gauge_field_thin[4];

// nHYP stuff
EXTERN Real alpha_smear[3];
EXTERN su3_matrix_f *hyplink1[4][4]; /* Needed for other stuff, too */
EXTERN su3_matrix_f *hyplink2[4][4];
EXTERN su3_matrix_f *Sigma[4];
EXTERN su3_matrix_f *SigmaH[4];
EXTERN su3_matrix_f *SigmaH2[4][4];

EXTERN su3_matrix_f *Staple1[4][4];
EXTERN su3_matrix_f *Staple2[4][4];
EXTERN su3_matrix_f *Staple3[4];

EXTERN su3_matrix_f *LambdaU[4];
EXTERN su3_matrix_f *Lambda1[4];
EXTERN su3_matrix_f *Lambda2[4];

/* renamed from tempmat1, for force_nhyp
   also used to compute gauge action */
EXTERN su3_matrix_f *tempmat_nhyp1;
EXTERN su3_matrix_f *tempmat_nhyp2;

#ifdef NHYP_JACOBI
EXTERN Matrix Qj, Vj;
#define JACOBI_HIST_MAX 10
EXTERN int jacobi_hist[JACOBI_HIST_MAX], jacobi_total;
EXTERN Real jacobi_avrg;
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
