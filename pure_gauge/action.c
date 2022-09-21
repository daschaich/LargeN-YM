// -----------------------------------------------------------------
// Measure total action, as needed by the hybrid Monte Carlo algorithm
#include "pg_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Helper function: Magnitude squared of an anti-hermitian matrix
Real ahmat_mag_sq(anti_hermitmat *pt) {
  register int i;
  register Real x, sum = 0.0;

  for (i = 0; i < NCOL; i++) {
    x = pt->im_diag[i];
    sum += 0.5 * x * x;
  }
  for (i = 0; i < N_OFFDIAG; i++) {
    x = pt->m[i].real;
    sum += x * x;
    x = pt->m[i].imag;
    sum += x * x;
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
  g_act = -beta * volume * (ssplaq + stplaq);
  h_act = hmom_action();
  tot = a * g_act + h_act;    // Fixed a=1 without LLR
  node0_printf("ACTION: g, h, tot = %8g %.8g %.8g\n", g_act, h_act, tot);
  return tot;
}
// -----------------------------------------------------------------
