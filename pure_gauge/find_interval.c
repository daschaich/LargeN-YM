// -----------------------------------------------------------------
// Find initial configuration in target energy interval
#include "pg_includes.h"
//#define FINDING_DEBUG

void findEint(double E_min) {
#ifdef LLR
  int Efound = -99;
  int counter = 0, count_max = 2000;
  double E, dtime, beta_sav = beta;
  double E_max = E_min + delta;
  node0_printf("Searching for energy interval [%.8g, %.8g]\n", E_min, E_max);

  dtime = -dclock();
  E = gauge_action();
  if (E >= E_min && E <= E_max)
    Efound = 1;

  while (Efound < 0 && counter < count_max) {
    update();
    E = gauge_action() * beta_sav / beta;
    if (E >= E_min && E <= E_max) {
      Efound = 1;
      beta = beta_sav;    // Reset beta
    }
    else if (E > E_max)
      beta += 0.1; // Increases plaq --> decreases E
    else  // E < Emin
      beta -= 0.1;
#ifdef FINDING_DEBUG
    node0_printf("FINDING E %.8g %.8g %d\n", E, beta, counter);
#endif
    node0_printf("FINDING E %.8g %.8g %d\n", E, beta, counter);
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
