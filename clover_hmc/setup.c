// -----------------------------------------------------------------
// Dynamical nHYP Wilson-clover setup
#include "cl_dyn_includes.h"
#define IF_OK if (status==0)

// Each node has a params structure for passing simulation parameters
#include "params.h"
params par_buf;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// On node zero, read lattice size, seed, nflavors and send to others
int initial_set() {
  int prompt, status;
  if (mynode() == 0) {
    // Print banner
    printf("SU(%d) with clover fermions, DIMF = %d, fermion rep = ",
           NCOL, DIMF);
#if FREP == fundamental
    printf("fundamental\n");
#elif FREP == symmetric2
    printf("two-index symmetric\n");
#elif FREP == antisymmetric2
    printf("two-index antisymmetric\n");
#else
    printf("unrecognized... shutting down\n");
    terminate(1);
#endif
    printf("Microcanonical simulation with refreshing\n");
    printf("Machine = %s, with %d nodes\n", machine_type(), numnodes());
#if SMEAR_LEVEL == 3
    printf("nHYP links, reading alpha_smear parameters from infile\n");
#elif SMEAR_LEVEL < 3
    printf("Smearing with %d levels (last %d unused)\n",
           SMEAR_LEVEL, 3 - SMEAR_LEVEL);
#endif
    printf("  IR_STAB = %.4g\n", (Real)IR_STAB);
#if NCOL == 3
    printf("  EPS_SQ = %.4g\n", (Real)EPS_SQ);
#elif NCOL == 4
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
  #if NCOL == 4
    printf("  TOL_NHYP = %.4g for SU(4)\n", (Real)TOL_NHYP);
  #endif
#endif
    printf("BETA_FREP turned on\n");
#ifdef HMC_ALGORITHM
    printf("Hybrid Monte Carlo algorithm\n");
#endif
#ifdef PHI_ALGORITHM
    printf("Phi algorithm\n");
#else
    printf("R algorithm\n");
#endif
    printf("LU preconditioning\n");
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
    IF_OK status += get_i(stdin, prompt, "nx", &par_buf.nx);
    IF_OK status += get_i(stdin, prompt, "ny", &par_buf.ny);
    IF_OK status += get_i(stdin, prompt, "nz", &par_buf.nz);
    IF_OK status += get_i(stdin, prompt, "nt", &par_buf.nt);
    IF_OK status += get_i(stdin, prompt, "iseed", &par_buf.iseed);

    if (status > 0)
      par_buf.stopflag = 1;
    else
      par_buf.stopflag = 0;
  }

  // Broadcast parameter buffer from node 0 to all other nodes
  broadcast_bytes((char *)&par_buf, sizeof(par_buf));
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
  one_ov_N = 1.0 / (Real)NCOL;
  volume = nx * ny * nz * nt;
  one_ov_vol = 1.0 / (Real)volume;
  total_iters = 0;
  return prompt;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Allocate all space for fields
void make_fields() {
  Real size = (Real)(8.0 * sizeof(su3_matrix_f));
  FIELD_ALLOC_VEC(gauge_field, su3_matrix_f, 4);
  FIELD_ALLOC_VEC(gauge_field_thin, su3_matrix_f, 4);

  /* move here alloc for clov? */
  size += (Real)(sizeof(su3_matrix));
  FIELD_ALLOC(f_mn, su3_matrix);

  // CG stuff
  size += (Real)((5.0 + 2.0 * MAX_MASSES) * sizeof(wilson_vector));
  FIELD_ALLOC(g_rand, wilson_vector);
  FIELD_ALLOC_VEC(chi, wilson_vector, MAX_MASSES);
  FIELD_ALLOC_VEC(psi, wilson_vector, MAX_MASSES);
  FIELD_ALLOC(p, wilson_vector);
  FIELD_ALLOC(mp, wilson_vector);
  FIELD_ALLOC(r, wilson_vector);
  FIELD_ALLOC(tempwvec, wilson_vector);
#ifdef PHI_ALGORITHM
  size += (Real)(MAX_MASSES * sizeof(wilson_vector));
  FIELD_ALLOC_VEC(old_psi, wilson_vector, MAX_MASSES);
#endif

  size += (Real)(8.0 * sizeof(half_wilson_vector));
  FIELD_ALLOC_VEC(htmp, half_wilson_vector, 8);

  // nHYP stuff
  size += (Real)(16.0 * sizeof(su3_matrix_f));
  FIELD_ALLOC_VEC(Sigma, su3_matrix_f, 4);
  FIELD_ALLOC_VEC(SigmaH, su3_matrix_f, 4);
  FIELD_ALLOC_VEC(Staple3, su3_matrix_f, 4);
  FIELD_ALLOC_VEC(LambdaU, su3_matrix_f, 4);

#if SMEAR_LEVEL > 1
  size += (Real)(28.0 * sizeof(su3_matrix_f));
  FIELD_ALLOC_VEC(Lambda1, su3_matrix_f, 4);
  FIELD_ALLOC_MAT_OFFDIAG(hyplink2, su3_matrix_f, 4);
  FIELD_ALLOC_MAT_OFFDIAG(Staple2, su3_matrix_f, 4);
#endif

#if SMEAR_LEVEL == 3
  size += (Real)(44.0 * sizeof(su3_matrix_f));
  FIELD_ALLOC_VEC(Lambda2, su3_matrix_f, 4);
  FIELD_ALLOC_MAT_OFFDIAG(hyplink1, su3_matrix_f, 4);
  FIELD_ALLOC_MAT_OFFDIAG(Staple1, su3_matrix_f, 4);
  FIELD_ALLOC_MAT(SigmaH2, su3_matrix_f, 4, 4);
#endif

  size += (Real)(3.0 * sizeof(su3_matrix_f));
  FIELD_ALLOC(tempmatf, su3_matrix_f);
  FIELD_ALLOC(tempmatf2, su3_matrix_f);
  FIELD_ALLOC(staplef, su3_matrix_f);

  size += (Real)(3.0 * sizeof(su3_matrix));
  FIELD_ALLOC(tempmat, su3_matrix);
  FIELD_ALLOC(tempmat2, su3_matrix);
  FIELD_ALLOC(staple, su3_matrix);

  size *= (Real)sites_on_node;
  node0_printf("Mallocing %.1f MBytes per core for fields\n", size / 1e6);

#if NCOL == 4
#ifdef NHYP_JACOBI
  Qj = AllocateMatrix(NCOL);
  Vj = AllocateMatrix(NCOL);
#endif
#endif
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int setup()   {
  int prompt;

  // Print banner, get volume, nflavors, seed
  prompt = initial_set();
  // Initialize the node random number generator
  initialize_prn(&node_prn, iseed, volume + mynode());
  // Initialize the layout functions, which decide where sites live
  setup_layout();
  // Allocate space for lattice, set up coordinate fields
  make_lattice();
  // Set up neighbor pointers and comlink structures
  make_nn_gathers();
  // Allocate space for fields
  make_fields();

  return prompt;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Read in parameters for SU(N) Monte Carlo
// prompt=1 indicates prompts are to be given for input
int readin(int prompt) {
  int status, i;
  Real x;

  // On node zero, read parameters and send to all other nodes
  if (this_node == 0) {
    printf("\n\n");
    status = 0;

    // Warms, trajecs
    IF_OK status += get_i(stdin, prompt, "warms", &par_buf.warms);
    IF_OK status += get_i(stdin, prompt, "trajecs", &par_buf.trajecs);
    IF_OK status += get_f(stdin, prompt, "traj_length", &par_buf.traj_length);

    // Number of pseudofermions
    IF_OK status += get_i(stdin, prompt, "number_of_PF", &par_buf.num_masses);
    if (par_buf.num_masses > MAX_MASSES || par_buf.num_masses < 1) {
      printf("num_masses = %d must be <= %d and >0!\n",
             par_buf.num_masses, MAX_MASSES);
      status++;
    }
    if (par_buf.num_masses != 2)
      printf("WARNING: code not tested for number_of_PF !=2\n");

    // Number of fermion steps
    for (i = 0; i < MAX_MASSES; i++)
      IF_OK status += get_i(stdin, prompt, "nstep", &par_buf.nsteps[i]);

    // Number of gauge steps
    IF_OK status += get_i(stdin, prompt, "nstep_gauge",
                          &par_buf.nsteps[MAX_MASSES]);

    // Trajectories between propagator measurements
    IF_OK status += get_i(stdin, prompt, "traj_between_meas",
                          &par_buf.propinterval);

    // beta, beta_frep, kappa and shift
    IF_OK status += get_f(stdin, prompt, "beta", &par_buf.beta);
    IF_OK status += get_f(stdin, prompt, "beta_frep", &par_buf.beta_frep);
    IF_OK status += get_f(stdin, prompt, "kappa", &par_buf.kappa);
    if (par_buf.num_masses > 1)
      IF_OK status += get_f(stdin, prompt, "shift", &par_buf.shift);

    // Clover coefficient and u0
    IF_OK status += get_f(stdin, prompt, "clov_c", &par_buf.clov_c);
    IF_OK status += get_f(stdin, prompt, "u0", &par_buf.u0);

    // Smearing parameters
    IF_OK status += get_f(stdin, prompt, "alpha_hyp0", &par_buf.alpha_hyp0);
    IF_OK status += get_f(stdin, prompt, "alpha_hyp1", &par_buf.alpha_hyp1);
    IF_OK status += get_f(stdin, prompt, "alpha_hyp2", &par_buf.alpha_hyp2);

    // Maximum numbers of conjugate gradient iterations and restarts
    IF_OK status += get_i(stdin, prompt, "max_cg_iterations", &par_buf.niter);
    IF_OK status += get_i(stdin, prompt, "max_cg_restarts", &par_buf.nrestart);

    // Error per site for update CG
    IF_OK {
      status += get_f(stdin, prompt, "error_per_site", &x);
      par_buf.rsqmin = x * x;
    }

    // Error per site for propagator measurement CG
    IF_OK {
      status += get_f(stdin, prompt, "error_for_propagator", &x);
      par_buf.rsqprop = x * x;
    }

    // Find out what kind of starting lattice to use
    IF_OK status += ask_starting_lattice(stdin,  prompt, &par_buf.startflag,
                                         par_buf.startfile);

    // Find out what to do with lattice at end
    IF_OK status += ask_ending_lattice(stdin,  prompt, &(par_buf.saveflag),
                                       par_buf.savefile);
    IF_OK status += ask_ildg_LFN(stdin,  prompt, par_buf.saveflag,
                                 par_buf.stringLFN);

    if (status > 0)
      par_buf.stopflag = 1;
    else
      par_buf.stopflag = 0;
  }

  // Broadcast parameter buffer from node0 to all other nodes
  broadcast_bytes((char *)&par_buf, sizeof(par_buf));
  if (par_buf.stopflag != 0)
    normal_exit(0);

  warms = par_buf.warms;
  trajecs = par_buf.trajecs;
  traj_length=par_buf.traj_length;

  num_masses = par_buf.num_masses;
  for (i = 0; i < MAX_MASSES + 1; i++)
    nsteps[i] = par_buf.nsteps[i];

  propinterval = par_buf.propinterval;

  beta = par_buf.beta;
  beta_frep = par_buf.beta_frep;
  kappa = par_buf.kappa;
  mkappa = -kappa;
  mkappaSq = -kappa * kappa;
  if (num_masses > 1)
    shift = par_buf.shift;
  else
    shift = 0.0;
  ishift = cmplx(0.0, shift);
  CNEGATE(ishift, mishift);

  clov_c = par_buf.clov_c;
  u0 = par_buf.u0;
  CKU0 = kappa * clov_c / (u0 * u0 * u0);

  alpha_smear[0] = par_buf.alpha_hyp0;
  alpha_smear[1] = par_buf.alpha_hyp1;
  alpha_smear[2] = par_buf.alpha_hyp2;

  niter = par_buf.niter;
  nrestart = par_buf.nrestart;
  rsqmin = par_buf.rsqmin;
  rsqprop = par_buf.rsqprop;

  startflag = par_buf.startflag;
#ifdef SPECTRUM
  fixflag = COULOMB_GAUGE_FIX;
#else
  fixflag = NO_GAUGE_FIX;
#endif
  saveflag = par_buf.saveflag;
  strcpy(startfile, par_buf.startfile);
  strcpy(savefile, par_buf.savefile);

  // Do whatever is needed to get lattice
  startlat_p = reload_lattice(startflag, startfile);

  // Create DIMFxDIMF matrices 'link' from NCOLxNCOL matrices 'linkf'
  fermion_rep();
  return 0;
}
// -----------------------------------------------------------------
