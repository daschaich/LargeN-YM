// -----------------------------------------------------------------
// Measure total action, as needed by the hybrid Monte Carlo algorithm
// When this routine is called the conjugate gradient should already
// have been run on the even sites, so that
// the vectors psi[i] contain (M^dag.M)^(-1) chi[i]
#include "cl_dyn_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Fermion contribution to the action
double fermion_action() {
  register int i;
  register site *s;
  double sum = 0.0;
  double shiftSq = shift * shift;

  FOREVENSITES(i, s) {
    if (num_masses == 1)
      wvec_rdot_sum(&(psi[0][i]), &(chi[0][i]), &sum);
    else {
      // Level 0 is the slow invert
      // Level 1 is the fast (shifted) invert
      wvec_rdot_sum(&(chi[0][i]), &(chi[0][i]), &sum);
      sum += shiftSq * (double)wvec_rdot(&(psi[0][i]), &(chi[0][i]));
      wvec_rdot_sum(&(psi[1][i]), &(chi[1][i]), &sum);
    }
  }
  g_doublesum(&sum);
  return sum;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Sum magnitude squared of an antihermition matrix
void ahmat_mag_sq(anti_hermitmat *pt, Real *c) {
  register Real x;

  x = pt->m00im;    *c += 0.5 * x * x;
  x = pt->m11im;    *c += 0.5 * x * x;
  x = pt->m01.real; *c += x * x;
  x = pt->m01.imag; *c += x * x;
#if NCOL > 2
  x = pt->m22im;    *c += 0.5 * x * x;
  x = pt->m02.real; *c += x * x;
  x = pt->m02.imag; *c += x * x;
  x = pt->m12.real; *c += x * x;
  x = pt->m12.imag; *c += x * x;
#if NCOL > 3
  x = pt->m33im;    *c += 0.5 * x * x;
  x = pt->m03.real; *c += x * x;
  x = pt->m03.imag; *c += x * x;
  x = pt->m13.real; *c += x * x;
  x = pt->m13.imag; *c += x * x;
  x = pt->m23.real; *c += x * x;
  x = pt->m23.imag; *c += x * x;
#endif
#endif
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Gauge momentum contribution to the action
double hmom_action() {
  register int i, dir;
  register site *s;
  double sum = 0.0;

  FORALLUPDIR(dir) {
    FORALLSITES(i, s)
      ahmat_mag_sq(&(s->mom[dir]), &sum);
  }
  g_doublesum(&sum);
  return sum;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
double action() {
  double g_act, g_frep = 0.0, h_act, f_act, tot;
  double ssplaq, stplaq, ssplaq_frep, stplaq_frep;

  plaquette(&ssplaq, &stplaq);
  ssplaq = 1.0 - ssplaq / (Real)NCOL;
  stplaq = 1.0 - stplaq / (Real)NCOL;
  // Three space--space and three space--time plaquette orientations
  g_act = beta * 3.0 * volume * (ssplaq + stplaq);

  if (fabs(beta_frep) > IMAG_TOL) {
    plaquette_frep(&ssplaq_frep, &stplaq_frep);
    ssplaq_frep = 1.0 - ssplaq_frep / (Real)DIMF;
    stplaq_frep = 1.0 - stplaq_frep / (Real)DIMF;
    g_frep = beta_frep * 3.0 * volume * (ssplaq_frep + stplaq_frep);
  }

  h_act = hmom_action();
  f_act = fermion_action();
  tot = g_act + g_frep + h_act + f_act;
  node0_printf("ACTION: g, rep, h, f, tot = %.8g %.8g %.8g %.8g %.8g\n",
               g_act, g_frep, h_act, f_act, tot);

  return tot;
}
// -----------------------------------------------------------------
