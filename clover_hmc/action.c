// -----------------------------------------------------------------
// Measure total action, as needed by the hybrid Monte Carlo algorithm
// When this routine is called the conjugate gradient should already
// have been run on the even sites, so that
// the vector psi contains (M^dag.M)^(-1) chi
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
      sum += (double)wvec_rdot(&(s->psi[0]), &(s->chi[0]));
    else{
      // Level 0 is the slow invert
      // Level 1 is the fast (shifted) invert
      sum += (double)wvec_rdot(&(s->chi[0]), &(s->chi[0]));
      sum += shiftSq * (double)wvec_rdot(&(s->psi[0]), &(s->chi[0]));
      sum += (double)wvec_rdot(&(s->psi[1]), &(s->chi[1]));
    }
  }
  g_doublesum(&sum);
  return sum;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Magnitude squared of an antihermition matrix
Real ahmat_mag_sq(anti_hermitmat *pt) {
  register Real x, sum;

  x = pt->m00im;    sum  = 0.5 * x * x;
  x = pt->m11im;    sum += 0.5 * x * x;
  x = pt->m01.real; sum += x * x;
  x = pt->m01.imag; sum += x * x;
#if (NCOL > 2)
  x = pt->m22im;    sum += 0.5 * x * x;
  x = pt->m02.real; sum += x * x;
  x = pt->m02.imag; sum += x * x;
  x = pt->m12.real; sum += x * x;
  x = pt->m12.imag; sum += x * x;
#if (NCOL > 3)
  x = pt->m33im;    sum += 0.5 * x * x;
  x = pt->m03.real; sum += x * x;
  x = pt->m03.imag; sum += x * x;
  x = pt->m13.real; sum += x * x;
  x = pt->m13.imag; sum += x * x;
  x = pt->m23.real; sum += x * x;
  x = pt->m23.imag; sum += x * x;
#endif
#endif
  return sum;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Gauge momentum contribution to the action
double hmom_action() {
  register int i,dir;
  register site *s;
  double sum = 0.0;

  FORALLUPDIR(dir) {
    FORALLSITES(i, s)
      sum += (double)ahmat_mag_sq(&(s->mom[dir]));
  }
  g_doublesum(&sum);
  return sum;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
double action() {
  double ssplaq, g_act = 0, g_frep = 0, h_act = 0, f_act = 0, tot;
#ifndef IMP
  double stplaq, ssplaq_frep, stplaq_frep;
#endif

#ifdef IMP
  gauge_action(&ssplaq);
  g_act = beta * ssplaq / (Real)NCOL;
#else
  plaquette(&ssplaq, &stplaq);
  ssplaq = 1.0 - ssplaq / (Real)NCOL;
  stplaq = 1.0 - stplaq / (Real)NCOL;
  // Three space--space and three space--time plaquette orientations
  g_act = beta * 3.0 * volume * (ssplaq + stplaq);

  if (fabs(beta_frep) > 1e-6) {
    plaquette_frep(&ssplaq_frep, &stplaq_frep);
    ssplaq_frep = 1.0 - ssplaq_frep / (Real)DIMF;
    stplaq_frep = 1.0 - stplaq_frep / (Real)DIMF;
    g_frep = beta_frep * 3.0 * volume * (ssplaq_frep + stplaq_frep);
  }
#endif

  h_act = hmom_action();
  f_act = fermion_action();
  tot = g_act + g_frep + h_act + f_act;
  node0_printf("ACTION: g, rep, h, f, tot = %.8g %.8g %.8g %.8g %.8g\n",
               g_act, g_frep, h_act, f_act, tot);

  return tot;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Copy NCOLxNCOL and DIMFxDIMF gauge fields
// as four-component arrays of su3_matrix_f and su3_matrix, respectively
void gauge_field_copy_f(field_offset src, field_offset dest) {
  register int i, dir, src2, dest2;
  register site *s;

  FORALLSITES(i, s) {
    src2 = src;
    dest2 = dest;
    FORALLUPDIR(dir) {
      su3mat_copy_f((su3_matrix_f *)F_PT(s, src2),
                    (su3_matrix_f *)F_PT(s, dest2));
      src2 += sizeof(su3_matrix_f);
      dest2 += sizeof(su3_matrix_f);
    }
  }
}

void gauge_field_copy(field_offset src, field_offset dest) {
  register int i, dir,src2,dest2;
  register site *s;

  FORALLSITES(i, s) {
    src2 = src;
    dest2 = dest;
    FORALLUPDIR(dir) {
      su3mat_copy((su3_matrix *)F_PT(s, src2),
                  (su3_matrix *)F_PT(s, dest2));
      src2 += sizeof(su3_matrix);
      dest2 += sizeof(su3_matrix);
    }
  }
}
// -----------------------------------------------------------------
