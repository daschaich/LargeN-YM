// -----------------------------------------------------------------
/* Version for dynamical clover fermions with Symanzik/tadople
  improved gauge field  */

#include "cl_dyn_includes.h"
#include <string.h>
int initial_set();
void make_fields();
#define IF_OK if (status==0)

/* Each node has a params structure for passing simulation parameters */
#include "params.h"
params par_buf;

int setup()   {
  int prompt;

  /* print banner, get volume, nflavors, seed */
  prompt=initial_set();
  /* initialize the node random number generator */
  initialize_prn(&node_prn,iseed,volume+mynode());
  /* Initialize the layout functions, which decide where sites live */
  setup_layout();
  /* allocate space for lattice, set up coordinate fields */
  make_lattice();
  /* set up neighbor pointers and comlink structures */
  make_nn_gathers();
  /* allocate space for fields */
  make_fields();

  return(prompt);
}


/* SETUP ROUTINES */
int initial_set() {
  int prompt,status;
  /* On node zero, read lattice size, seed, nflavors and send to others */
  if (mynode() == 0) {
    /* print banner */
/* stringification kludge from gnu preprocessor manual
   http://gcc.gnu.org/onlinedocs/cpp/Stringification.html */
#define XSTR(s) STR(s)
#define STR(s) #s
/* end kludge */
    printf("SU(%d) with clover fermions, DIMF = %d, fermion rep = "
           XSTR(FREP) "\n",NCOL,DIMF);
#ifdef BETA_FREP
//  printf("BETA_FREP turned on\n");
#endif
#ifdef NHYP
    printf("nHYP links with %d smearing level(s)\n", SMEAR_LEVEL);
    printf("  Reading alpha_smear parameters from infile\n");
#if (SMEAR_LEVEL < 3)
    printf("  REMEMBER: last %d smearing parameter(s) not used in this run\n",
           3 - SMEAR_LEVEL);
#endif
    printf("  IR_STAB = %e\n", (Real)IR_STAB);
#if (NCOL == 3)
    printf("  EPS_SQ = %e\n", (Real)EPS_SQ);
#elif (NCOL == 4)
  #ifdef NHYP_JACOBI
    printf("Using Jacobi\n  TOL_JACOBI = %e\n  MAX_JACOBI_ITERS = %d\n",
     TOL_JACOBI, MAX_JACOBI_ITERS);
  #else
    printf("  EPS_SQ_4 = %e\n  EPS_SQ_3 = %e\n",
            (Real)EPS_SQ_4, (Real)EPS_SQ_3);
  #endif
#endif
#ifdef NHYP_DEBUG
    printf("NHYP_DEBUG turned on\n");
#if (NCOL == 4)
    printf("  TOL_NHYP = %e\n", TOL_NHYP);
  #ifdef NHYP_JACOBI
    printf("  PRINT_JACOBI_ITERS = %d\n", PRINT_JACOBI_ITERS);
  #else
    printf("  TOL_ACOS = %.4g\n", (Real)TOL_ACOS);
  #endif
#endif
#endif
#endif /* NHYP */
    printf("\nMicrocanonical simulation with refreshing\n");
    printf("MIMD version 7ish\n");
    printf("Machine = %s, with %d nodes\n",machine_type(),numnodes());
#ifdef HMC_ALGORITHM
    printf("Hybrid Monte Carlo algorithm\n");
#endif
#ifdef PHI_ALGORITHM
    printf("Phi algorithm\n");
#else
    printf("R algorithm\n");
#endif
#ifdef LU
    printf("LU preconditioning\n");
#endif
#ifdef SPECTRUM
    printf("With spectrum measurements\n");
#endif
    time_stamp("start");

    status = get_prompt(stdin,  &prompt);

    IF_OK status += get_i(stdin, prompt, "nflavors", &par_buf.nflavors);
#ifdef PHI_ALGORITHM
    IF_OK if (par_buf.nflavors != 2) {
      printf("Error: Use phi algorithm only for two flavors\n");
      terminate(-1);
    }
#endif

    IF_OK status += get_i(stdin, prompt,"nx", &par_buf.nx);
    IF_OK status += get_i(stdin, prompt,"ny", &par_buf.ny);
    IF_OK status += get_i(stdin, prompt,"nz", &par_buf.nz);
    IF_OK status += get_i(stdin, prompt,"nt", &par_buf.nt);

    IF_OK status += get_i(stdin, prompt,"iseed", &par_buf.iseed);

    if (status>0) par_buf.stopflag=1; else par_buf.stopflag=0;
  } /* end if (mynode()==0) */

  /* Node 0 broadcasts parameter buffer to all other nodes */
  broadcast_bytes((char *)&par_buf,sizeof(par_buf));

  if (par_buf.stopflag != 0)
    normal_exit(0);

  nflavors = par_buf.nflavors;
  nx = par_buf.nx;
  ny = par_buf.ny;
  nz = par_buf.nz;
  nt = par_buf.nt;
  iseed = par_buf.iseed;

  this_node = mynode();
  number_of_nodes = numnodes();
  volume = nx * ny * nz * nt;
  total_iters = 0;
  return prompt;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Read in parameters for SU(N) Monte Carlo
int readin(int prompt) {
  // prompt=1 indicates prompts are to be given for input
  int status, i;
  Real x;

  // On node zero, read parameters and send to all other nodes
  if (this_node==0) {
    printf("\n\n");
    status = 0;

    // Warms, trajecs
    IF_OK status += get_i(stdin, prompt,"warms", &par_buf.warms);
    IF_OK status += get_i(stdin, prompt,"trajecs", &par_buf.trajecs);
    IF_OK status += get_f(stdin, prompt,"traj_length", &par_buf.traj_length);

       // Number of pseudofermions
        IF_OK status += get_i(stdin, prompt,"number_of_PF", &par_buf.num_masses);
        if (par_buf.num_masses>MAX_MASSES || par_buf.num_masses<1) {
            printf("num_masses = %d must be <= %d and >0!\n", par_buf.num_masses, MAX_MASSES);
            status++;
        }
        if (par_buf.num_masses != 2)printf("WARNING: code not tested for number_of_PF !=2\n");
        for (i=0;i<MAX_MASSES;i++)
    IF_OK status += get_i(stdin, prompt, "nstep", &par_buf.nsteps[i]);
        IF_OK status += get_i(stdin, prompt, "nstep_gauge", &par_buf.nsteps[MAX_MASSES]);
    /* trajectories between propagator measurements */
    IF_OK status += get_i(stdin, prompt,"traj_between_meas",
        &par_buf.propinterval);



    /* get couplings and broadcast to nodes */
    /* beta, kappa */
    IF_OK status += get_f(stdin, prompt,"beta", &par_buf.beta);
#ifdef BETA_FREP
    IF_OK status += get_f(stdin, prompt,"beta_frep", &par_buf.beta_frep);
#endif
    IF_OK status += get_f(stdin, prompt,"kappa", &par_buf.kappa);
    if (par_buf.num_masses > 1)
    IF_OK status += get_f(stdin, prompt,"shift", &par_buf.shift);

    // Clover coefficient and u0
    IF_OK status += get_f(stdin, prompt, "clov_c", &par_buf.clov_c);
    IF_OK status += get_f(stdin, prompt, "u0", &par_buf.u0);

#ifdef NHYP
    // Smearing parameters
    IF_OK status += get_f(stdin, prompt, "alpha_hyp0", &par_buf.alpha_hyp0);
    IF_OK status += get_f(stdin, prompt, "alpha_hyp1", &par_buf.alpha_hyp1);
    IF_OK status += get_f(stdin, prompt, "alpha_hyp2", &par_buf.alpha_hyp2);
#endif

    // Maximum numbers of conjugate gradient iterations and restarts
    IF_OK status += get_i(stdin, prompt, "max_cg_iterations", &par_buf.niter);
    IF_OK status += get_i(stdin, prompt, "max_cg_restarts", &par_buf.nrestart);

    // Error per site for update CG
    IF_OK {
      status += get_f(stdin, prompt,"error_per_site", &x);
      par_buf.rsqmin = x * x;
    }

    // Error per site for propagator measurement CG
    IF_OK {
      status += get_f(stdin, prompt,"error_for_propagator", &x);
      par_buf.rsqprop = x * x;
    }

    /* find out what kind of starting lattice to use */
    IF_OK status += ask_starting_lattice(stdin,  prompt, &par_buf.startflag,
  par_buf.startfile);

    /* find out what to do with lattice at end */
    IF_OK status += ask_ending_lattice(stdin,  prompt, &(par_buf.saveflag),
           par_buf.savefile);
    IF_OK status += ask_ildg_LFN(stdin,  prompt, par_buf.saveflag,
          par_buf.stringLFN);

    if (status > 0)par_buf.stopflag=1; else par_buf.stopflag=0;
  } /* end if (this_node==0) */

  broadcast_bytes((char *)&par_buf,sizeof(par_buf));
  if (par_buf.stopflag != 0)
    normal_exit(0);

  warms = par_buf.warms;
  trajecs = par_buf.trajecs;
  traj_length=par_buf.traj_length;
  num_masses = par_buf.num_masses;
  for (i = 0; i < MAX_MASSES + 1; i++)
    nsteps[i] = par_buf.nsteps[i];
  if (num_masses > 1)
    shift = par_buf.shift;

  propinterval = par_buf.propinterval;
  startflag = par_buf.startflag;
  saveflag = par_buf.saveflag;
  niter = par_buf.niter;
  nrestart = par_buf.nrestart;
  rsqmin = par_buf.rsqmin;
  rsqprop = par_buf.rsqprop;

  beta = par_buf.beta;
#ifdef BETA_FREP
  beta_frep = par_buf.beta_frep;
#endif
  kappa = par_buf.kappa;
  clov_c = par_buf.clov_c;
  u0 = par_buf.u0;
#ifdef NHYP
  alpha_smear[0] = par_buf.alpha_hyp0;
  alpha_smear[1] = par_buf.alpha_hyp1;
  alpha_smear[2] = par_buf.alpha_hyp2;
#endif
  strcpy(startfile,par_buf.startfile);
  strcpy(savefile,par_buf.savefile);

  // Load part of inversion control structure for generic inverters
  // Note different convention: linear resid rather than quadratic rsq
  qic.min = 0;
  qic.max = niter;
  qic.nrestart = nrestart;
  qic.resid = sqrt(rsqprop);

  // Load part of Dirac matrix parameters
  dcp.Kappa = kappa;
  dcp.Clov_c = clov_c;
  dcp.U0 = u0;

  // Do whatever is needed to get lattice
  startlat_p = reload_lattice(startflag, startfile);
  /* generates the fermion's irrep lattice,
     including spatially twisted b.c. for SF */
  fermion_rep();
  return 0;
}

/* allocate all space for fields */
void make_fields() {

int memfield;

/* move here alloc for clov? */

#ifdef NHYP
    FIELD_ALLOC_VEC(gauge_field,su3_matrix_f,4);
    FIELD_ALLOC_VEC(gauge_field_thin,su3_matrix_f,4);

    FIELD_ALLOC_VEC(Sigma,su3_matrix_f,4);
    FIELD_ALLOC_VEC(SigmaH,su3_matrix_f,4);
    FIELD_ALLOC_VEC(Staple3,su3_matrix_f,4);
    FIELD_ALLOC_VEC(LambdaU,su3_matrix_f,4);
    memfield=26;
#if (SMEAR_LEVEL>1)
    FIELD_ALLOC_VEC(Lambda1,su3_matrix_f,4);
    FIELD_ALLOC_MAT_OFFDIAG(hyplink2,su3_matrix_f,4);
    FIELD_ALLOC_MAT_OFFDIAG(Staple2,su3_matrix_f,4);
    memfield=54;
#endif

#if (SMEAR_LEVEL==3)
    FIELD_ALLOC_VEC(Lambda2,su3_matrix_f,4);
    FIELD_ALLOC_MAT_OFFDIAG(hyplink1,su3_matrix_f,4);
    FIELD_ALLOC_MAT_OFFDIAG(Staple1,su3_matrix_f,4);
    FIELD_ALLOC_MAT(SigmaH2,su3_matrix_f,4,4);
    memfield=98;
#endif

    FIELD_ALLOC(tempmat_nhyp1,su3_matrix_f);
    FIELD_ALLOC(tempmat_nhyp2,su3_matrix_f);

    node0_printf("Mallocing %.1f MBytes per node for fields\n",
            (double)sites_on_node * memfield * sizeof(su3_matrix_f)/1e6);

#if (NCOL==4)
#ifdef NHYP_JACOBI
    Qj = AllocateMatrix(NCOL);
    Vj = AllocateMatrix(NCOL);
#endif
#endif

#endif
}
// -----------------------------------------------------------------
