// -----------------------------------------------------------------
// Measure total action, as needed by the hybrid Monte Carlo algorithm
#include "pg_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Helper function: Magnitude squared of an anti-hermitian matrix
Real ahmat_mag_sq(anti_hermitmat *ah) {
  register int i;
  register Real sum;

  sum = ah->im_diag[0] * ah->im_diag[0];
  for (i = 1; i < NCOL; i++)
    sum += ah->im_diag[i] * ah->im_diag[i];
  sum *= 0.5;

  for (i = 0; i < N_OFFDIAG; i++) {
    sum += ah->m[i].real * ah->m[i].real;
    sum += ah->m[i].imag * ah->m[i].imag;
  }

  return sum;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Gauge momentum contribution to the action
double hmom_action() {
#ifdef HMC
  register int i, dir;
  register site *s;
  double sum = 0.0;

  FORALLSITES(i,s) {
    FORALLUPDIR(dir)
      sum += (double)ahmat_mag_sq(&(s->mom[dir]));
  }
  g_doublesum(&sum);
  return sum;
#else
  return -99.0;
#endif
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
double gauge_action() {
  double ssplaq, stplaq;

  plaquette(&ssplaq, &stplaq);
  ssplaq = ssplaq * one_ov_N;
  stplaq = stplaq * one_ov_N;
  // Three space--space and three space--time plaquette orientations
  return (-beta * 3.0 * volume * (ssplaq + stplaq));
}

double action(double E_min) {
  double g_act, h_act, tot;

  g_act = gauge_action();
  node0_printf("ACTION: g, h, ");
#ifdef LLR
  double td = 0.0, w_act = 0.0;
  if (constrained == 1) {
    // Add gaussian window contribution
    //td = g_act - E_min - 0.5 * delta;
    //w_act = exp(-0.5 * td * td / deltaSq);
    //w_act = 0.5 * td * td / deltaSq;
  }

#endif

  h_act = hmom_action();
  node0_printf("tot = %.8g %.8g ", g_act, h_act);

  tot = g_act + h_act;
#ifdef LLR
  if (constrained == 1)
    tot = a * g_act + h_act;
//    tot = a * g_act + h_act + w_act;
#endif
  node0_printf("%.8g\n", tot);
  return tot;
}
// -----------------------------------------------------------------
