#ifndef _LATTICE_H
#define _LATTICE_H
/****************************** lattice.h ********************************/

/* include file for MIMD QCD program, version 2
   This file defines global scalars and the fields in the lattice. */

#include "defines.h"
#include "../include/generic_quark_types.h"
#include "../include/macros.h"    /* For MAXFILENAME */
#include "../include/io_lat.h"    /* For gauge_file */
#include "../include/random.h"    /* For double_prn */

/* Begin definition of site structure */

#include "../include/su3.h"
#include "../include/random.h"   /* For double_prn */

#ifdef NHYP_JACOBI
#include "../include/jacobi.h"
#endif

/* not in use */
#ifdef NHYP_GMP
#include "../include/gmp_stuff.h"
#endif

/* The lattice is an array of sites.  */
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
 	su3_matrix_f old_linkf[4];
	/* For accept/reject */
#endif

	/* antihermitian momentum matrices in each direction */
 	anti_hermitmat mom[4];

	/* Wilson complex vectors */
	wilson_vector g_rand;	/* gaussian random vector */
	wilson_vector psi[MAX_MASSES];	/* solution vector */
        wilson_vector chi[MAX_MASSES];	/* source vector */
#ifdef BI
        wilson_vector sss;      /* internal biconjugate gradient vector */
        wilson_vector ttt;      /* internal biconjugate gradient vector */
        wilson_vector rv;       /* internal biconjugate gradient vector */
        wilson_vector v;        /* internal biconjugate gradient vector */
#endif
        wilson_vector p;        /* conjugate gradient change vector */
        wilson_vector mp;       /* another CG vector */
        wilson_vector tmp;      /* another temporary CG vector */
	wilson_vector r; 	/* residue */
	wilson_vector invp;     /* used in paca_[tx] */
	wilson_vector inva;     /* used in paca_[tx] */

        /* wilson half vector (temporary used in dslash_w_site) */
        half_wilson_vector htmp[8];
/**half_wilson_vector htmp2[8];**/ /* TEMP FOR TESTING */
#ifdef PHI_ALGORITHM
 	wilson_vector old_psi[MAX_MASSES];	/* For predicting next psi */
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
  /* a bit wasteful of space but needed for Hasenbusch, I think */
  su3_matrix Force[4];


    /* align to double word boundary */
	/**double space2;**/
} site;

/* End definition of site structure */

/* Definition of globals */

#ifdef CONTROL
#define EXTERN
#else
#define EXTERN extern
#endif

/* The following are global scalars */
EXTERN	int nx,ny,nz,nt;	/* lattice dimensions */
EXTERN  int volume;		/* volume of lattice = nx*ny*nz*nt */
EXTERN	int iseed;		/* random number seed */
EXTERN	int warms,trajecs,niter,nrestart,propinterval,nflavors;
EXTERN    float traj_length;

/*  masses (ultimately for hasenbusch...  */
EXTERN  int num_masses;         /* max number of masses <= MAX_MASSES */
EXTERN Real shift;
EXTERN    int nsteps[MAX_MASSES+1];



EXTERN	Real rsqmin,rsqprop,beta,kappa,clov_c,u0;
#ifdef BETA_FREP
EXTERN	Real beta_frep;
#endif
EXTERN	Real epsilon;
EXTERN	char startfile[MAXFILENAME],savefile[MAXFILENAME];
EXTERN  double g_ssplaq, g_stplaq;
EXTERN  double_complex linktrsum;
EXTERN  u_int32type nersc_checksum;
EXTERN  char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable **/
EXTERN	int startflag;	/* beginning lattice: CONTINUE, RELOAD, FRESH */
EXTERN  int fixflag;  /* gauge fix: COULOMB_GAUGE_FIX, NO_GAUGE_FIX */
EXTERN	int saveflag;	/* do with lattice: 1=save; */
EXTERN	int total_iters;
/* boundary-flip switch, made consistent for higher-rep code */
EXTERN	int current_boundary;
EXTERN	int current_boundary_x;

/* Some of these global variables are node dependent */
/* They are set in "make_lattice()"                  */
EXTERN	int sites_on_node;		/* number of sites on this node */
EXTERN	int even_sites_on_node;	/* number of even sites on this node */
EXTERN	int odd_sites_on_node;	/* number of odd sites on this node */
EXTERN	int number_of_nodes;	/* number of nodes in use */
EXTERN  int this_node;		/* node number of this node */

EXTERN  int debugflag;

/* Schrodinger Functional stuff  */
#ifdef SF
/* this is to allow keeping the original flag in the code,
   to distinguish original MILC from our modifications                */
#define SCHROED_FUN

/* additional globals            */
/* spatially constant boundary fields and their eta-derivatives
   NOTE: spatial boundary links at t==0 stored also in *lattice       */
EXTERN  su3_matrix   link_bndr_up[3],   link_driv_up[3],
                     link_bndr_dn[3],   link_driv_dn[3];
EXTERN  su3_matrix_f linkf_bndr_up[3],  linkf_driv_up[3],
                     linkf_bndr_dn[3],  linkf_driv_dn[3];
#ifdef NHYP
EXTERN  su3_matrix_f linkf_zero[3];
#endif
/* fermion phase factors for twisted b.c. in x,y,z directions         */
EXTERN  complex ferm_phases[3];

/* gauge action:
   perturbative weight for space-time plaquettes at boundaries        */
EXTERN  Real c_t;
EXTERN  int bc_flag;          /* flag for gauge field boundary values */
EXTERN  int ct_flag;          /* flag for activating perturbative c_t */
EXTERN  Real ferm_twist_phase; /* twist angle for spatial b.c.        */

/* SF flag options.  Must be defined even when not doing SF  */
/* bc_flag options */
#define TRIVIAL   0
#define ALPHA     1  /* for SU(2) only */
#define ALPHA_ONE 2  /* for SU(3) only */
#define ALPHA_TWO 3  /* for SU(3) only */
#define SF_SYM    4  /* now for SU(4) only */
/* ct_flag options */
#define TREE_LEVEL 0
#define ONE_LOOP   1 /* for SU(3) with fundamental irrep only */
#define TWO_LOOP   2 /* for SU(3) with fundamental irrep only */

#endif /* SF */

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

/* definitions for arccos */
/* nobody seems to use this
#define NMAX 401

EXTERN Real acos_table[NMAX];
EXTERN Real acos_deriv[NMAX];
EXTERN Real delta_acos;
*/

EXTERN quark_invert_control qic;
EXTERN dirac_clover_param dcp;

EXTERN gauge_file *startlat_p;
EXTERN gauge_file *savelat_p;

/* The lattice is a single global variable - (actually this is the
   part of the lattice on this node) */
EXTERN site *lattice;

/* Vectors for addressing */
/* Generic pointers, for gather routines */
#define N_POINTERS 8	/* Number of generic pointers */
/* NEED 8 WHEN GAUGEFIXING; 10 in original nhyp code */
EXTERN char ** gen_pt[N_POINTERS];

#ifdef NHYP
/* hard-coding of alpha_smear controlled by defines.h */
#ifndef HARD_CODE_SMEAR
EXTERN Real alpha_smear[3];
#endif

/* control for re-use of fermion_force in sf_coupling                */
/* during update, set it to FORCE                                    */
/* set it to SF_COUPLING when calculating sf_coupling with fat links */
EXTERN int sf_coupling_flag;
#define FORCE 747
#define SF_COUPLING 767

/* global scalar, accummulates the fat-link contribution to sf_coupling */
EXTERN double f_fat;

/* prototype fields for NHYP links */
EXTERN su3_matrix_f *gauge_field[4];
EXTERN su3_matrix_f *gauge_field_thin[4];

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
/* EXTERN su3_matrix_f *save_gf[4]; */

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

/* gmp not in use */
/* allocate variables globally for gmp (cf include/gmp_stuff) */
#ifdef NHYP_GMP
/* variables needed by helper routines (generic_nhyp/gmp_helpers.c) */
EXTERN mpf_t m_tmp, m_tr;
EXTERN msu3_matrix_f m_Q, m_Qtmp1, m_Qtmp2;
/* variables used in nhyp_SU4.c */
EXTERN mpf_t m_c0,m_c1,m_c2,m_c3,m_b0,m_b1,m_b2,m_S,m_R,m_c0sq;
#endif

#endif /* NHYP */

/* timing */
/* up to 20 concurrent timers */
#ifdef TIMING
EXTERN double tmptime[20];
/* congrad_cl_m from: coupling, grsource, update_o, use tmptime[0] */
EXTERN double time_dcongrad;
/* fermion_force from: coupling, update_h_sf2, update_o, use tmptime[1] */
EXTERN double time_fermion_force;
/* fermion_rep from: setup, update_o, use tmptime[2] */
EXTERN double time_fermion_rep;
/* block_nhyp from: fermion_rep, use tmptime[3] */
EXTERN double time_block_nhyp;
/* compute_fhb from: block_nhyp, force_nhyp, use tmptime[4] */
EXTERN double time_compute_fhb;
/* use tmptime[5] */
EXTERN double time_gmp;
/* use tmptime[6] */
EXTERN double time_jacobi;
#endif

#endif /* _LATTICE_H */
