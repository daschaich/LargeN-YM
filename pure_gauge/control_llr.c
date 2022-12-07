// -----------------------------------------------------------------
// Main procedure and helpers for pure-gauge over-relaxed quasi-heatbath
#define CONTROL
#include "pg_includes.h"

int main(int argc, char *argv[]) {
  int prompt;
  int traj_done, RMcount;//, Nmeas = 0;
  double ss_plaq, st_plaq, E, dtime, save_a;
  double Reweightexpect;    // Reweighted expectation value of the energy

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

#ifndef HMC
  // Monitor overall acceptance in monteconst_e.c
  accept = 0;
  reject = 0;
  double rate;
#endif

  // Check: compute initial plaquette and energy
  plaquette(&ss_plaq, &st_plaq);
  E = gauge_action();
  node0_printf("START %.8g %.8g %.8g %.8g\n",
               ss_plaq, st_plaq, ss_plaq + st_plaq, E);

  // Unconstrained warmup sweeps before searching for energy interval
  for (traj_done = 0; traj_done < warms; traj_done++)
    update();
  node0_printf("WARMUPS COMPLETED\n");

  // Terminates if interval not found
  constrained = 0;
  findEint();

  // Robbins--Monro (RM) iterations
  constrained = 1;
  for (RMcount = 0; RMcount < ait; RMcount++) {
    // Constrained warm-up sweeps in each RM iteration, with a=1
    save_a = a;
    a = 1.0;
    for (traj_done = 0; traj_done < warms; traj_done++)
      updateconst_e();
    a = save_a;

    Reweightexpect = 0.0;
    for (traj_done = 0; traj_done < trajecs; traj_done++) {
      updateconst_e();
      // Accumulate after update
      Reweightexpect += gauge_action();

      // More expensive measurements every "measinterval" sweeps
      if ((traj_done % measinterval) == (measinterval - 1)) {
//        Nmeas++;
        // Nothing yet...
      }
    }
    Reweightexpect /= trajecs;
    Reweightexpect -= Emin + 0.5 * delta;

    // Hard-code under-relaxation to begin after 100 RM iterations
    if (RMcount < 100)
      a += 12.0 * Reweightexpect / deltaSq;
    else
      a += 12.0 * Reweightexpect / (deltaSq * (RMcount - 99));
    node0_printf("RM ITER %d a %.8g\n", RMcount + 1, a);
    // TODO: I think acceptance rate for each RM iteration
    //       would be more interesting than the overall one below...
  }
  node0_printf("RUNNING COMPLETED\n");

  // Check: compute final plaquette and energy
  plaquette(&ss_plaq, &st_plaq);
  E = gauge_action();
  node0_printf("STOP %.8g %.8g %.8g %.8g\n",
               ss_plaq, st_plaq, ss_plaq + st_plaq, E);

#ifndef HMC
  rate = (double)accept / ((double)(accept + reject));
  node0_printf("Overall acceptance %d of %d = %.4g\n",
               accept, accept + reject, rate);
#endif
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
