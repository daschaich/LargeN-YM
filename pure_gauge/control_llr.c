// -----------------------------------------------------------------
// Main procedure and helpers for pure-gauge over-relaxed quasi-heatbath
#define CONTROL
#include "pg_includes.h"

int main(int argc, char *argv[]) {
  int prompt;
  int traj_done, RMcount, Ncount, Nint, Intcount;
  double ss_plaq, st_plaq, E, dtime;
  double RestrictEV;    // Restricted expectation value <<E - E_i>>

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

  // Set up number of intervals
  if (Emax < Emin) {
    node0_printf("ERROR: Emax smaller than Emin\n");
    terminate(1);
  }
  Nint = (int)((Emax - Emin) / delta) + 1;
  node0_printf("Nint %d\n", Nint);
  Real Eint[Nint], aint[Nint];

#ifndef HMC
  // Monitor quasi-heatbath acceptance in monteconst_e.c
  accept = 0;
  reject = 0;
#endif

  // Cycle over intervals in this job
  for (Intcount = 0; Intcount < Nint; Intcount++) {
    aint[Intcount] = 0;
    Eint[Intcount] = Emin + Intcount * delta;

    // Cycle over jackknife samples in this interval
    for (Ncount = 0; Ncount < Nj; Ncount++) {
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
      // Otherwise sets initial guess for a
      constrained = 0;
      findEint(Eint[Intcount]);

      // Newton--Raphson (NR) and Robbins--Monro (RM) iterations
      constrained = 1;
      for (RMcount = 0; RMcount < (NRiter + RMiter); RMcount++) {
        // Constrained warm-up sweeps in each iteration
        for (traj_done = 0; traj_done < warms; traj_done++)
          updateconst_e(Eint[Intcount]);

        RestrictEV = 0.0;
        for (traj_done = 0; traj_done < trajecs; traj_done++) {
          updateconst_e(Eint[Intcount]);
          // Accumulate after update
          RestrictEV += gauge_action();
        }
        RestrictEV /= trajecs;
        RestrictEV -= Eint[Intcount] + 0.5 * delta;

        // Implement initial NR iterations
        // Replace deltaSq/12 --> deltaSq to control fluctuations
        if (RMcount < NRiter) {
          if (fabs(RestrictEV) < a_cut)
            a += RestrictEV / deltaSq;
          else if (RestrictEV > 0.0)
            a += a_cut / deltaSq;
          else
            a -= a_cut / deltaSq;
        }
        else {    // Otherwise under-relax
          if (fabs(RestrictEV) < a_cut)
            a += RestrictEV / (deltaSq *  (RMcount - NRiter + 1));
          else if (RestrictEV > 0.0)
            a += a_cut / (deltaSq * (RMcount - NRiter + 1));
          else
            a -= a_cut / (deltaSq * (RMcount - NRiter + 1));
        }
        node0_printf("RM ITER %d RestrictEV %.8g a %.8g\n",
                     RMcount + 1, RestrictEV, a);
      }
      aint[Intcount] += a;

      // Reload lattice for next jackknife sample
      if (Ncount < Nj - 1) {
        node0_printf("Resetting lattice\n");
        startlat_p = reload_lattice(startflag, startfile);
      }
    }
    aint[Intcount] /= (double)Nj;
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
