// -----------------------------------------------------------------
// Main procedure for Wilson flow
#define CONTROL
#include "wflow_includes.h"

int main(int argc, char *argv[])  {
  int prompt;
  double dtime;

  // Set up
  setlinebuf(stdout); // DEBUG
  initialize_machine(&argc, &argv);
  g_sync();
  if (remap_stdio_from_args(argc, argv) == 1)
    terminate(1);

  prompt = setup();
  if (readin(prompt) != 0) {
    node0_printf("ERROR in readin, aborting\n");
    terminate(1);
  }
  dtime = -dclock();

  // Set up loop tables
  if (num_meas > 0)
    make_loop_table();

  // Run Wilson flow!
  // Adaptive step size stuff all handled in this routine
  wflow(S, A);

  node0_printf("RUNNING COMPLETED\n");
  dtime += dclock();
  node0_printf("Time = %.4g seconds\n", dtime);
  fflush(stdout);

  // Save lattice if requested
  if (saveflag != FORGET)
    save_lattice(saveflag, savefile, stringLFN);

  normal_exit(0);         // Needed by at least some clusters
  return 0;
}
// -----------------------------------------------------------------
