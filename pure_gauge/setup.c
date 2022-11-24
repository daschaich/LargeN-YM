// -----------------------------------------------------------------
// SU(N) pure-gauge setup
#include "pg_includes.h"
#define IF_OK if (status == 0)

// Each node has a params structure for passing simulation parameters
#include "params.h"
params par_buf;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// On node zero, read lattice size, seed, nflavors and send to others
int initial_set() {
  int prompt = 0, status = 0;
  if (mynode() == 0) {
    // Print banner
#ifndef LLR
#ifndef HMC
    printf("SU(%d) pure-gauge over-relaxed heatbath algorithm\n", NCOL);
#else
    printf("SU(%d) pure-gauge hybrid Monte Carlo algorithm\n", NCOL);
#endif
#else
    printf("SU(%d) pure-gauge log-linear relaxation algorithm ", NCOL);
#ifndef HMC
    printf("using over-relaxed heatbath updates\n");
#else
    printf("using hybrid Monte Carlo updates\n");
#endif
#endif
    printf("Machine = %s, with %d nodes\n", machine_type(), numnodes());
    time_stamp("start");

    status = get_prompt(stdin, &prompt);
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

  nx = par_buf.nx;
  ny = par_buf.ny;
  nz = par_buf.nz;
  nt = par_buf.nt;
  iseed = par_buf.iseed;

  this_node = mynode();
  number_of_nodes = numnodes();
  volume = nx * ny * nz * nt;
  one_ov_vol = 1.0 / (Real)volume;
  one_ov_N = 1.0 / (Real)NCOL;
  return prompt;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Allocate all space for fields
void make_fields() {
  Real size = (Real)(2.0 * sizeof(matrix));
  FIELD_ALLOC(tempmat, matrix);
  FIELD_ALLOC(tempmat2, matrix);

  size *= (Real)sites_on_node;
  node0_printf("Mallocing %.1f MBytes per core for fields\n", size / 1e6);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int setup() {
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
// Read in parameters for SU(N) pure-gauge evolution
// prompt=1 indicates prompts are to be given for input
int readin(int prompt) {
  int status;

  // On node zero, read parameters and send to all other nodes
  if (this_node == 0) {
    printf("\n\n");
    status = 0;

    // Warms, sweeps
    IF_OK status += get_i(stdin, prompt, "warms", &par_buf.warms);
    IF_OK status += get_i(stdin, prompt, "trajecs", &par_buf.trajecs);

    // Trajectories between more expensive measurements
    IF_OK status += get_i(stdin, prompt, "traj_between_meas",
                          &par_buf.measinterval);

    // beta
    IF_OK status += get_f(stdin, prompt, "beta", &par_buf.beta);

#ifndef HMC
    // Over-relaxation steps per sweep
    IF_OK status += get_i(stdin, prompt, "ora_steps", &par_buf.ora_steps);

    // Quasi-heatbath steps per sweep
    IF_OK status += get_i(stdin, prompt, "qhb_steps", &par_buf.qhb_steps);
#else
    // HMC steps per sweep
    IF_OK status += get_i(stdin, prompt, "hmc_steps", &par_buf.hmc_steps);

    // HMC trajectory length
    IF_OK status += get_f(stdin, prompt, "traj_length", &par_buf.traj_length);
#endif

#ifdef LLR
    // LLR stuff
    // Energy range to scan (min to max)
    IF_OK status += get_f(stdin, prompt, "Emin", &par_buf.Emin);

    // Size of energy interval delta
    IF_OK status += get_f(stdin, prompt, "delta", &par_buf.delta);

    // Starting value for LLR observable 'a'
    IF_OK status += get_f(stdin, prompt, "a", &par_buf.a);

    // Number of iterations for Robbins--Monro algorithm
    IF_OK status += get_i(stdin, prompt, "ait", &par_buf.ait);
#endif

    // Find out what kind of starting lattice to use
    IF_OK status += ask_starting_lattice(stdin, prompt, &par_buf.startflag,
                                         par_buf.startfile);

    // Find out what to do with lattice at end
    IF_OK status += ask_ending_lattice(stdin, prompt, &(par_buf.saveflag),
                                       par_buf.savefile);
    IF_OK status += ask_ildg_LFN(stdin, prompt, par_buf.saveflag,
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
  measinterval = par_buf.measinterval;
  beta = par_buf.beta;

#ifndef HMC
  ora_steps = par_buf.ora_steps;
  qhb_steps = par_buf.qhb_steps;
#else
  hmc_steps = par_buf.hmc_steps;
  traj_length = par_buf.traj_length;
#endif

#ifdef LLR
  a = par_buf.a;
  Emin = par_buf.Emin * volume;
  delta = par_buf.delta * volume;
  Emax = Emin + delta;
  deltaSq = delta * delta;
  ait = par_buf.ait;
#endif

  startflag = par_buf.startflag;
  saveflag = par_buf.saveflag;
  strcpy(startfile, par_buf.startfile);
  strcpy(savefile, par_buf.savefile);
  strcpy(stringLFN, par_buf.stringLFN);

  // Allocate some more arrays to be used by LAPACK
  // in generic/reunitarize.c
  Rwork = malloc(sizeof *Rwork * 5 * NCOL);
  eigs = malloc(sizeof *eigs * NCOL);
  store = malloc(sizeof *store * 2 * NCOL * NCOL);
  work = malloc(sizeof *work * 6 * NCOL);
  junk = malloc(sizeof *junk * NCOL);
  left = malloc(sizeof *left * 2 * NCOL * NCOL);
  right = malloc(sizeof *right * 2 * NCOL * NCOL);

  // Do whatever is needed to get lattice
  startlat_p = reload_lattice(startflag, startfile);
  return 0;
}
// -----------------------------------------------------------------
