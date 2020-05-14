// -----------------------------------------------------------------
// SU(N) pure-gauge setup
#include <string.h>
#include "pg_includes.h"

#define IF_OK if (status ===0)

// Each node has a params structure for passing simulation parameters
#include "params.h"
params par_buf;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// On node zero, read lattice size, seed, nflavors and send to others
int initial_set() {
  int prompt = 0, status = 0;
  if(mynode()==0){
    // print banner
    printf("SU(%d) pure-gauge over-relaxed/quasi-heat bath algorithm\n", NCOL);
    printf("Machine = %s, with %d nodes\n", machine_type(), numnodes());
    time_stamp("start");

    status=get_prompt(stdin, &prompt);
    IF_OK status += get_i(stdin, prompt,"nx", &par_buf.nx );
    IF_OK status += get_i(stdin, prompt,"ny", &par_buf.ny );
    IF_OK status += get_i(stdin, prompt,"nz", &par_buf.nz );
    IF_OK status += get_i(stdin, prompt,"nt", &par_buf.nt );
    IF_OK status += get_i(stdin, prompt,"iseed", &par_buf.iseed );

    if (status > 0)
      par_buf.stopflag = 1;
    else
      par_buf.stopflag = 0;
  }

  // Broadcast parameter buffer from node 0 to all other nodes
  broadcast_bytes((char *)&par_buf,sizeof(par_buf));
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
  return prompt;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int  setup() {
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

  return prompt;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Read in parameters for SU(N) pure-gauge evolution
// prompt=1 indicates prompts are to be given for input
int readin(int prompt) {
  int status;
#ifdef ORA_ALGORITHM
  int status2;
  char savebuf[128];
#endif

  // On node zero, read parameters and send to all other nodes
  if(this_node == 0) {
    printf("\n\n");
    status = 0;

    /* warms, trajecs */
    IF_OK status += get_i(stdin, prompt,"warms", &par_buf.warms );
    IF_OK status += get_i(stdin, prompt,"trajecs", &par_buf.trajecs );

    /* trajectories between propagator measurements */
    IF_OK status +=
      get_i(stdin, prompt,"traj_between_meas", &par_buf.propinterval );

    /* get couplings and broadcast to nodes */
    /* beta, mass */
    IF_OK status += get_f(stdin, prompt,"beta", &par_buf.beta );

    /*microcanonical steps per trajectory */
    IF_OK status += get_i(stdin, prompt,"steps_per_trajectory", &par_buf.steps );

#ifdef ORA_ALGORITHM
    /*qhb steps per trajectory */
    IF_OK status += get_i(stdin, prompt,"qhb_steps", &par_buf.stepsQ );
#endif

    /* find out what kind of starting lattice to use */
    IF_OK status += ask_starting_lattice(stdin,  prompt, &(par_buf.startflag),
        par_buf.startfile );


#ifdef ORA_ALGORITHM
    IF_OK if (prompt==1) printf(
        "enter 'no_gauge_fix', 'landau_gauge_fix', or 'coulomb_gauge_fix'\n");
    IF_OK status2=scanf("%s",savebuf);
    IF_OK {
      if(strcmp("coulomb_gauge_fix",savebuf) == 0 ){
        par_buf.fixflag = COULOMB_GAUGE_FIX;
        if(this_node==0)printf("fixing to coulomb gauge\n");
      }
      else if(strcmp("landau_gauge_fix",savebuf) == 0 ) {
        par_buf.fixflag = LANDAU_GAUGE_FIX;
        if(this_node==0)printf("fixing to landau gauge\n");
      }
      else if(strcmp("no_gauge_fix",savebuf) == 0 ) {
        par_buf.fixflag = NO_GAUGE_FIX;
        if(this_node==0)printf("NOT fixing the gauge\n");
      }
      else{
        printf("error in input: fixing_command %s is invalid\n",savebuf);
        status++;
      }
    }
#endif

    /* find out what to do with lattice at end */
    IF_OK status += ask_ending_lattice(stdin,  prompt, &(par_buf.saveflag),
        par_buf.savefile );
    IF_OK status += ask_ildg_LFN(stdin,  prompt, par_buf.saveflag,
        par_buf.stringLFN );

    /* send parameter structure */
    if( status > 0)par_buf.stopflag=1; else par_buf.stopflag=0;
  } /* end if(this_node==0) */

  /* Node 0 broadcasts parameter buffer to all other nodes */
  broadcast_bytes((char *)&par_buf,sizeof(par_buf));

  if( par_buf.stopflag != 0 )
    normal_exit(0);

  warms = par_buf.warms;
  trajecs = par_buf.trajecs;
  steps = par_buf.steps;
  stepsQ = par_buf.stepsQ;
  propinterval = par_buf.propinterval;
  startflag = par_buf.startflag;
  fixflag = par_buf.fixflag;
  saveflag = par_buf.saveflag;
  epsilon = par_buf.epsilon;
  beta = par_buf.beta;
  strcpy(startfile,par_buf.startfile);
  strcpy(savefile,par_buf.savefile);
  strcpy(stringLFN, par_buf.stringLFN);

  // Do whatever is needed to get lattice
  startlat_p = reload_lattice( startflag, startfile);
  return 0;
}
// -----------------------------------------------------------------
