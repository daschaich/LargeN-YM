// -----------------------------------------------------------------
// Main procedure and helpers for pure-gauge over-relaxed quasi-heatbath
#define CONTROL
#include "pg_includes.h"
//#define DEBUG_PRINT
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Find initial configuration in desired energy interval
void findEint(double Eint, double delta) {
  bool Efound = false;
  int counter = 0, count_max = 2000;
  double ss_plaq, st_plaq, dtime, energy, energyref, betaref = beta;

  node0_printf("Searching for energy interval [%.8g, %.8g]\n",
               Eint, Eint + delta);

  dtime = -dclock();
  energy = action(&ss_plaq, &st_plaq);
  if (energy >= Eint && energy <= (Eint + delta))
  {
    Efound = true;
    beta = betaref;
  }

  while(Efound == false && counter < count_max) {
    update();
    energyref = action(&ss_plaq, &st_plaq) * betaref / beta;
    if (energyref >= Eint && energyref <= (Eint + delta))
    {
      Efound = true;
      beta = betaref;
    }
    else if (energyref>(Eint+delta))
      beta += 0.1;
    else  // energyref < Eint
      beta -= 0.1;
#ifdef DEBUG_PRINt
    node0_printf("energyref %.8g %.8g %d\n", energyref, beta, counter);
#endif

    counter++;
  }

  dtime += dclock();
  if (Efound) {
    node0_printf("Energy %.4g in interval found ", energy);
    node0_printf("after %d updates in %.4g seconds\n", counter, dtime);
  }
  else {
    node0_printf("ERROR: Energy interval not found ");
    node0_printf("after %d updates in %.4g seconds\n", counter, dtime);
    node0_printf("Ended up with energy %.4g\n", energy);
    terminate(1);
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int main(int argc, char *argv[]) {
  bool Einterval = false;
  int prompt;
  int swp_done, Nint, Eint, jcount, acounter;
  int RMcount;      // Count Robbins--Monro iterations
  double ss_plaq, st_plaq, dtime, energy;

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

  // Check: compute initial plaquette and energy
  energy = action(&ss_plaq, &st_plaq);
  node0_printf("START %.8g %.8g %.8g %.8g\n",
               ss_plaq, st_plaq, ss_plaq + st_plaq, energy);

  srand(iseed); // !!!TODO: Shouldn't be needed...

  // Declarations that depend on input
  // Number of energy intervals, include interval starting at Emax
  Nint = (int)((Emax - Emin) / delta) + 1;
  double x0[Nint];  // Lower end of energy interval
  double a[Nint], a_i[Njacknife], a_i_new;
  double Reweightexpect; // Reweighted expectation value of the energy
  double *measurement = malloc(sweeps * sizeof(double));

  for(Eint=0;Eint<Nint;Eint++)
  {
    x0[Eint] = (double)(Eint)*delta + Emin;
    a[Eint] = 0.0;
    for(jcount = 0;jcount<Njacknife;jcount++)
    {
      a_i[jcount] = 1.0;
      a_i_new = 1.0;
      RMcount = 0;
      Einterval = false;

      // Perform warmup sweeps before searching for energy interval
      coldlat();
      for (swp_done = 0; swp_done < warms; swp_done++)
        update();
      node0_printf("%d WARMUPS COMPLETED FOR jcount=%d\n", warms, jcount);

      for(acounter=0;acounter<ait;acounter++)
      {
        a_i[jcount] = a_i_new;

        if (Einterval == false)
        {
          // Terminates if interval not found
          findEint(x0[Eint],delta);
          Einterval = true;

          for (swp_done = 0; swp_done < (2*warms); swp_done++)
            updateconst_e(x0[Eint], 1.0);
        }

        for (swp_done = 0; swp_done < warms; swp_done++)
          updateconst_e(x0[Eint], 1.0);

        // TODO: Compute variance from <O> and <O^2>
        Reweightexpect = 0;
        for (swp_done = 0; swp_done < sweeps; swp_done++) {
          measurement[swp_done] = action(&ss_plaq, &st_plaq);
          Reweightexpect += measurement[swp_done];
          updateconst_e(x0[Eint], a_i[jcount]);
        }
        Reweightexpect /= sweeps;

        Reweightexpect -= x0[Eint] + 0.5*delta;
        if (RMcount<100)
        {
          a_i_new = a_i[jcount] + 12/(delta*delta)*Reweightexpect;
          node0_printf("a = %.4g off\n", a_i_new);
        }
        else
        {
          a_i_new = a_i[jcount] + 12/(delta*delta*(RMcount+1-100))*Reweightexpect;
          node0_printf("a = %.4g\n", a_i_new);
        }
        RMcount++;
      }
      a[Eint] += a_i_new / Njacknife;
    }
    //node0_printf("a = %.4g \n", a[Eint]);
    //printf("%d \n",1.0);
  }
  node0_printf("RUNNING COMPLETED\n");

  // Check: compute final plaquette and energy
  energy = action(&ss_plaq, &st_plaq);
  node0_printf("STOP %.8g %.8g %.8g %.8g\n",
               ss_plaq, st_plaq, ss_plaq + st_plaq, energy);

  dtime += dclock();
  node0_printf("Time = %.4g seconds\n", dtime);
  fflush(stdout);

  // Save lattice if requested
  if (saveflag != FORGET)
    save_lattice(saveflag, savefile, stringLFN);

  free(measurement);
  normal_exit(0);         // Needed by at least some clusters
  return 0;
}
// -----------------------------------------------------------------
