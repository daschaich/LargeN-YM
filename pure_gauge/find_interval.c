// -----------------------------------------------------------------
// Find initial configuration in target energy interval
#include "pg_includes.h"
//#define FINDING_DEBUG

void findEint() {
#ifdef LLR
  int Efound = -99;
  int counter = 0, count_max = 2000;
  double E, dtime, beta_sav = beta;

  node0_printf("Searching for energy interval [%.8g, %.8g]\n", Emin, Emax);

  dtime = -dclock();
  E = gauge_action();
  if (E >= Emin && E <= Emax)
    Efound = 1;

  while (Efound < 0 && counter < count_max) {
    update();
    E = gauge_action() * beta_sav / beta;
    if (E >= Emin && E <= Emax) {
      Efound = 1;
      beta = beta_sav;    // Reset beta
    }
    else if (E > Emax)
      beta += 0.1;        // Increases plaq --> decreases E
    else  // E < Emin
      beta -= 0.1;
#ifdef FINDING_DEBUG
    node0_printf("FINDING E %.8g %.8g %d\n", E, beta, counter);
#endif

    counter++;
  }

  dtime += dclock();
  if (Efound > 0) {
    node0_printf("Energy %.8g in interval found ", E);
    node0_printf("after %d updates in %.4g seconds\n", counter, dtime);
  }
  else {
    node0_printf("ERROR: Energy interval not found ");
    node0_printf("after %d updates in %.4g seconds\n", counter, dtime);
    node0_printf("Ended up with energy %.4g\n", E);
    terminate(1);
  }
#endif
}
// -----------------------------------------------------------------
