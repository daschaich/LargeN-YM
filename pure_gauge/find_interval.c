// -----------------------------------------------------------------
// Find initial configuration in target energy interval
#include "pg_includes.h"

// Convert from plaquette to energy
double energy(double *ss_plaq, double *st_plaq) {
  double E;
  plaquette(ss_plaq, st_plaq);
  E = -3.0 * beta * volume * (*ss_plaq + *st_plaq);
  return E;
}

void findEint() {
#ifdef LLR
  int Efound = -99;
  int counter = 0, count_max = 2000;
  double E, ss_plaq, st_plaq, dtime, beta_sav = beta;

  node0_printf("Searching for energy interval [%.8g, %.8g]\n", Emin, Emax);

  dtime = -dclock();
  E = energy(&ss_plaq, &st_plaq);
  if (E >= Emin && E <= Emax)
    Efound = 1;

  while (Efound < 0 && counter < count_max) {
    update();
    E = energy(&ss_plaq, &st_plaq) * beta_sav / beta;
    if (E >= Emin && E <= Emax) {
      Efound = 1;
      beta = beta_sav;    // Reset beta
    }
    else if (E > Emax)
      beta += 0.1;
    else  // E < Emin
      beta -= 0.1;
#ifdef DEBUG_PRINt
    node0_printf("FINDING E %.8g %.8g %d\n", energy, beta, counter);
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
