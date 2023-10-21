// -----------------------------------------------------------------
// Main procedure and helpers for pure-gauge over-relaxed quasi-heatbath
#define CONTROL
#include "pg_includes.h"

int main(int argc, char *argv[]) {
  int prompt;
  int traj_done, RMcount, Ncount, nrintervals, Intcount;
  double ss_plaq, st_plaq, E, dtime, save_a, rate;
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

  if (Emax < Emin) {
    node0_printf("ERROR: Emax smaller than Emin\n");
    terminate(1);
  }
  nrintervals = (int)((Emax-Emin) / delta) + 1;
  node0_printf("nrintervals %d\n", nrintervals);
  double Eint[nrintervals];
  double aint[nrintervals];
  // Monitor overall acceptance in monteconst_e.c
  accept = 0;
  reject = 0;

  // Check: compute initial plaquette and energy
  plaquette(&ss_plaq, &st_plaq);
  E = gauge_action();
  node0_printf("START %.8g %.8g %.8g %.8g\n",
               ss_plaq, st_plaq, ss_plaq + st_plaq, E);
  save_a = a;
  for (Intcount = 0; Intcount < nrintervals; Intcount++) {
    aint[Intcount] = 0;
    Eint[Intcount] = Emin + Intcount*delta;
    for (Ncount = 0; Ncount < Njacknife; Ncount++) {
      a = save_a;
      constrained = 0;
      startlat_p = reload_lattice(startflag, startfile);
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
      a = 1.0;
      findEint(Eint[Intcount]);

      // Robbins--Monro (RM) iterations
      constrained = 1;
      for (RMcount = 0; RMcount < ait; RMcount++) {
        // Constrained warm-up sweeps in each RM iteration
        for (traj_done = 0; traj_done < warms; traj_done++)
          updateconst_e(Eint[Intcount]);

        Reweightexpect = 0.0;
        for (traj_done = 0; traj_done < trajecs; traj_done++) {
          updateconst_e(Eint[Intcount]);
          // Accumulate after update
          Reweightexpect += gauge_action();
        }
        Reweightexpect /= trajecs;
        Reweightexpect -= Eint[Intcount] + 0.5 * delta;

        // Hard-code under-relaxation to begin after initial NR iterations
        if (RMcount < NRiters) {
          if (fabs(Reweightexpect) < a_cut)
            a += Reweightexpect / deltaSq;
          else if (Reweightexpect > 0.0)
            a += a_cut / deltaSq;
          else
            a -= a_cut / deltaSq;
        }
        else {
          if (fabs(Reweightexpect) < a_cut)
            a += Reweightexpect / (deltaSq *  (RMcount - NRiters + 1));
          else if (Reweightexpect > 0.0)
            a += a_cut / (deltaSq * (RMcount - NRiters + 1));
          else
            a -= a_cut / (deltaSq * (RMcount - NRiters + 1));
        }
        node0_printf("RM ITER %d Reweightexpect %.8g a %.8g\n",
                     RMcount + 1, Reweightexpect, a);
      }
      aint[Intcount] = aint[Intcount] + a/((double)(Njacknife));
    }
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
