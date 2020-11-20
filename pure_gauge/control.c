// -----------------------------------------------------------------
// Main procedure for pure-gauge evolution
// Removed hybrid Monte Carlo updates,
// keeping just over-relaxed quasi-heat bath
#define CONTROL
#include "pg_includes.h"

int main(int argc, char *argv[]) {
  int traj_done;//, Nmeas = 0;
  int prompt;
  double ss_plaq, st_plaq, dtime;
  complex plp = cmplx(99.0, 99.0);

  // Set up
  setlinebuf(stdout); // DEBUG
  initialize_machine(&argc,&argv);
  g_sync();
  if (remap_stdio_from_args(argc, argv) == 1)
    terminate(1);

  // Load input and run (loop removed)
  prompt = setup();
  if (readin(prompt) != 0) {
    node0_printf("ERROR in readin, aborting\n");
    terminate(1);
  }
  dtime = -dclock();

  // Check: compute initial plaquette
  plaquette(&ss_plaq, &st_plaq);
  node0_printf("START %.8g %.8g %.8g\n", ss_plaq, st_plaq, ss_plaq + st_plaq);

  // Perform warmup trajectories
  for (traj_done = 0; traj_done < warms; traj_done++)
    update();
  node0_printf("WARMUPS COMPLETED\n");

  // Perform trajectories, reunitarizations and measurements
  for (traj_done = 0; traj_done < trajecs; traj_done++) {
    update();

    // Measure and print Polyakov loop and plaquette
    // after every trajectory
    plaquette(&ss_plaq, &st_plaq);
    plp = ploop(TUP);
    node0_printf("GMES %.8g %.8g %.8g %.8g\n",
                 plp.real, plp.imag, ss_plaq, st_plaq);
    fflush(stdout);

    // More expensive measurements every "propinterval" trajectories
    // Measure every "propinterval" trajectories
    if ((traj_done % propinterval) == (propinterval - 1)) {
//      Nmeas++;
      // Nothing yet...
    }
  }
  node0_printf("RUNNING COMPLETED\n");

  // Check: compute final plaquette
  plaquette(&ss_plaq, &st_plaq);
  node0_printf("STOP %.8g %.8g %.8g\n", ss_plaq, st_plaq, ss_plaq + st_plaq);

  dtime += dclock();
  node0_printf("Time = %.4g seconds\n", dtime);
  fflush(stdout);

  // Save lattice if requested
  if (saveflag != FORGET)
    save_lattice(saveflag, savefile, stringLFN);

  normal_exit(0);
  g_sync();         // Needed by at least some clusters
  return 0;
}
// -----------------------------------------------------------------
