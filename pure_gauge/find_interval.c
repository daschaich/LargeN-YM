// -----------------------------------------------------------------
// Find initial configuration in target energy interval
#include "pg_includes.h"

void findEint(double E_min) {
#ifdef LLR
  int Efound = -99, counter = 0;
  Real tr, dtime, beta_sav = beta;
  double E, E_max = E_min + delta;
  node0_printf("Searching for energy interval [%.8g, %.8g]\n", E_min, E_max);

  dtime = -dclock();
  E = gauge_action();
  node0_printf("FINDING E %.8g %.8g %d\n", E, beta, counter);
  if (E >= E_min && E <= E_max)
    Efound = 1;

  while (Efound < 0 && counter < FIND_MAX) {
    update_ora();
    E = gauge_action() * beta_sav / beta;
    if (E >= E_min && E <= E_max) {
      Efound = 1;
      a = beta;
      beta = beta_sav;    // Reset beta
    }
    else if (E > E_max) {
      counter++;
      tr = fabs(E - E_max);
      if (tr > 1.0)
        beta += 0.1;
      else if (tr > 0.1)
        beta += 0.01;
      else if (tr > 0.01)
        beta += 0.001;
      else
        beta += 0.0001;
    }
    else {                // E < Emin
      counter++;
      tr = fabs(E - E_min);
      if (tr > 1.0)
        beta -= 0.1;
      else if (tr > 0.1)
        beta -= 0.01;
      else if (tr > 0.01)
        beta -= 0.001;
      else
        beta -= 0.0001;
    }
#ifdef DEBUG_PRINt
    node0_printf("FINDING E %.8g %.8g %d\n", E, beta, counter);
#else
    if (counter % 10 == 0)
      node0_printf("FINDING E %.8g %.8g %d\n", E, beta, counter);
#endif
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
