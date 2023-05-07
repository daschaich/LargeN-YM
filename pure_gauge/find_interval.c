// -----------------------------------------------------------------
// Find initial configuration in target energy interval
#include "pg_includes.h"

void findEint(double E_min) {
#ifdef LLR
  int Efound = -99;
  int counter = 0, count_max = 2000;
  double E, dtime, beta_sav = beta;
  double E_max = E_min + delta;
  node0_printf("Searching for energy interval [%.8g, %.8g]\n", E_min, E_max);
#ifdef HMC
  int hmc_steps_save = hmc_steps;
  hmc_steps = 60;
  
  ora_steps = 5;
  qhb_steps = 5;
#endif
  
  dtime = -dclock();
  E = gauge_action();
  if (E >= E_min && E <= E_max)
    Efound = 1;

  while (Efound < 0 && counter < count_max) {
#ifndef HMC
    update();
#endif
#ifdef HMC
    relax();
    monte();
#endif
    E = gauge_action() * beta_sav / beta;
    if (E >= E_min && E <= E_max) {
      Efound = 1;
      a=beta;
      beta = beta_sav;    // Reset beta
    }
    else if (E > E_max){
      if(abs(E-(E_max))>1.0)
      {
        beta += 0.1;
      }
      else if (abs(E-E_max)>0.1)
      {
        beta += 0.01;
      }
      else if (abs(E-E_max)>0.01)
      {
        beta += 0.001;
      }
      else
      {
        beta += 0.0001;
      }
    }
    else{  // E < Emin
      if(abs(E-(E_min))>1.0)
      {
        beta -= 0.1;
      }
      else if (abs(E-E_min)>0.1)
      {
        beta -= 0.01;
      }
      else if (abs(E-E_min)>0.01)
      {
        beta -= 0.001;
      }
      else
      {
        beta -= 0.0001;
      }
    }
#ifdef DEBUG_PRINt
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
#ifdef HMC
  hmc_steps = hmc_steps_save;
#endif
#endif
}
// -----------------------------------------------------------------
