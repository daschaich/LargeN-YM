// -----------------------------------------------------------------
// Main procedure for pure-gauge over-relaxed quasi-heatbath
#define CONTROL
#include "pg_includes.h"

// -----------------------------------------------------------------
// Find initial configuration in desired energy interval
#ifdef HMCLLR
void findEint(double Eint, double delta) {
  bool Efound = false;
  int counter = 0, count_max = 2000;
  double ss_plaq, st_plaq, dtime, energy, energyref, betaref = beta;

  node0_printf("Searching for energy interval [%.8g, %.8g]\n",
               Eint, Eint + delta);

  dtime = -dclock();
  energy = U_action();
  if (energy >= Eint && energy <= (Eint + delta))
  {
    Efound = true;
    beta = betaref;
  }

  while(Efound == false && counter < count_max) {
    update();
    energyref = U_action() * betaref / beta;
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
    node0_printf("Energy %.4g in interval found ", energyref);
    node0_printf("after %d updates in %.4g seconds\n", counter, dtime);
  }
  else {
    node0_printf("ERROR: Energy interval not found ");
    node0_printf("after %d updates in %.4g seconds\n", counter, dtime);
    node0_printf("Ended up with energy %.4g\n", energyref);
    terminate(1);
  }
}
#endif
// -----------------------------------------------------------------

int main(int argc, char *argv[]) {
  int swp_done;//, Nmeas = 0;
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
  coldlat();
  // Check: compute initial plaquette
  plaquette(&ss_plaq, &st_plaq);
  node0_printf("START %.8g %.8g %.8g\n", ss_plaq, st_plaq, ss_plaq + st_plaq);
#ifdef HMCLLR
  findEint(Emin,delta);
#endif
  // Perform warmup sweeps
  for (swp_done = 0; swp_done < warms; swp_done++) {
#ifdef HMCLLR
    updateconst_e(Emin,1.0);
#else
    update();
#endif
   } 
  node0_printf("WARMUPS COMPLETED\n");

  // Perform sweeps, reunitarizations and measurements
  for (swp_done = 0; swp_done < sweeps; swp_done++) {
#ifdef HMCLLR
    updateconst_e(Emin,1.0);
#else
    update();
#endif
    // Measure and print Polyakov loop and plaquette after every sweep
    plaquette(&ss_plaq, &st_plaq);
    plp = ploop(TUP);
    node0_printf("GMES %.8g %.8g %.8g %.8g\n",
                 plp.real, plp.imag, ss_plaq, st_plaq);
    fflush(stdout);

    // More expensive measurements every "measinterval" sweeps
    if ((swp_done % measinterval) == (measinterval - 1)) {
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

  normal_exit(0);         // Needed by at least some clusters
  return 0;
}
// -----------------------------------------------------------------
