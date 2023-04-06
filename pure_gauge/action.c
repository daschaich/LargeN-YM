// -----------------------------------------------------------------
// Measure total action, as needed by the hybrid Monte Carlo algorithm
#include "pg_includes.h"
//#define DEBUG_PRINT
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
#ifdef DEBUG_PRINT
  node0_printf("Computing gauge_action\n");
#endif

  plaquette(&ssplaq, &stplaq);
  ssplaq = 1.0 - ssplaq * one_ov_N;
  stplaq = 1.0 - stplaq * one_ov_N;
#if 0 // TESTING
  ssplaq *= one_ov_N;
  stplaq *= one_ov_N;
#endif
  // Three space--space and three space--time plaquette orientations
  return (beta * 3.0 * volume * (ssplaq + stplaq));
}

double action(double E_min) {
  double g_act, h_act, tot;

  g_act = gauge_action();
  h_act = hmom_action();
#ifdef LLR
  g_act *= a;

  // TODO: Best to add window action here,
  //       then accept/reject step can stay standardized
//  if (constrained == 1)
//    tot = a * g_act + h_act + w_act;
#endif
  tot = g_act + h_act;
  node0_printf("ACTION: g, h, tot = %.8g %.8g %.8g\n", g_act, h_act, tot);
  return tot;
}
// -----------------------------------------------------------------
