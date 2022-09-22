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
  register int i, dir;
  register site *s;
  double sum = 0.0;

  FORALLSITES(i,s) {
    FORALLUPDIR(dir)
      sum += (double)ahmat_mag_sq(&(s->mom[dir]));
  }
  g_doublesum(&sum);
  return sum;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
double action() {
  double ssplaq, stplaq, g_act, h_act, tot;

  plaquette(&ssplaq, &stplaq);
  ssplaq = 1.0 - ssplaq * one_ov_N;
  stplaq = 1.0 - stplaq * one_ov_N;
  // Three space--space and three space--time plaquette orientations
  g_act = beta * 3.0 * volume * (ssplaq + stplaq);

  h_act = hmom_action();
  tot = g_act + h_act;
  //tot = a * g_act + h_act;    // TODO: Fix a=1 without LLR
  node0_printf("ACTION: g, h, tot = %8g %.8g %.8g\n", g_act, h_act, tot);
  return tot;
}
// -----------------------------------------------------------------
