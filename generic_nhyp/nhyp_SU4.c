// -----------------------------------------------------------------
// Helper routines for SU(4) nHYP blocking

#if NCOL != 4
  #error "Wrong NCOL in nhyp.c!"
#endif

#include "Power.h"
#include "jacobi.c"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute the Cayley Hamilton coefficients for the inverse square root
//   1 / sqrt(Q) = f[0] I + f[1] Q +f[2] Q^2 + f[3] Q^3
// The b[][] matrix contains the derivatives of f wrt to the traces of Q
//   b[i][j] = d f[i] / d c[j] with c[j] = 1/(j+1) trace Q^(j+1)
// The flag compute_b switches on/off the computation of these derivatives
#ifndef NHYP_DEBUG
void compute_fhb(matrix_f *Q, Real *f, Real b[NCOL][NCOL], int compute_b)
#else
void compute_fhb(matrix_f *Omega, matrix_f *Q,
                 Real *f, Real b[NCOL][NCOL], int compute_b)
#endif

{
  double sg0, sg1, sg2, sg3;
  double u, v, w, x, den;
  int i, j;

#ifdef TIMING
  TIC(4)
#endif

  /**** Find eigenvalues of Q *******************************************/
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      Qj.M[i][j].real=Q->e[i][j].real;
      Qj.M[i][j].imag=Q->e[i][j].imag;
    }
  }

  Jacobi(&Qj, &Vj, TOL_JACOBI);

  sg0 = Qj.M[0][0].real;
  sg1 = Qj.M[1][1].real;
  sg2 = Qj.M[2][2].real;
  sg3 = Qj.M[3][3].real;

#ifdef NHYP_DEBUG
  /* all eigenvalues must be positive */
  if (sg0 < 0 || sg1 < 0 || sg2 < 0 || sg3 < 0) {
    printf("NHYP_DEBUG_FATAL: negative eigenvalue\n sg0 = %.12e\n sg1 = %.12e\n sg2 = %.12e\n sg3 = %.12e\n",
        sg0, sg1, sg2, sg3);
    terminate(1);
  }
#endif

  /* square roots of the eigenvalues */
  sg0 = sqrt(sg0);
  sg1 = sqrt(sg1);
  sg2 = sqrt(sg2);
  sg3 = sqrt(sg3);

  /* The coefficients and derivatives are symmetric under permuations
     of g0, g1, g2. Since they are also polynomial in the sqrt's, they can be
     reduced to elementary symmetric functions u, v, w, x
     */
  u = sg0+sg1+sg2+sg3;
  v = sg0*(sg1+sg2+sg3)+sg1*(sg2+sg3)+sg2*sg3;
  w = sg0*sg1*(sg2+sg3)+(sg0+sg1)*sg2*sg3;
  x = sg0*sg1*sg2*sg3;

  /* calculate f's and their derivatives */
  den = -(x*(-(u*v*w) + w*w + u*u*x));
  f[0] = -w*w*w - u*u*w*x + v*w*x + u*(v*w*w - v*v*x + x*x);
  f[0] = f[0] / den;
  f[1] = 2*u*u*v*w - u*(v*v*v + 2*w*w) - u*u*u*x + w*(v*v + x);
  f[1] = f[1] / den;
  f[2] = u*u*u*v - u*u*w + 2*v*w + u*(-2*v*v + x);
  f[2] = f[2] / den;
  f[3] = -(u*v) + w;
  f[3] = f[3] / den;

  /* most of the time, the coefficients are all we need */
  if (compute_b == 0) {
#ifdef TIMING
    TOC(4, time_compute_fhb)
#endif
    return;
  }

  den =  2*Power(den, 3);
  /* Cf.: den = -2*x*x*x*Power(-(u*v*w) + w*w + u*u*x, 3);  */

  b[0][0] = Power(u, 7)*Power(x, 5)
    + Power(u, 6)*w*Power(x, 3)*(Power(w, 2) - 8*v*x)
    + Power(u, 5)*Power(x, 2)*(-3*v*Power(w, 4)
        + 16*Power(v, 2)*Power(w, 2)*x - 4*Power(v, 3)*Power(x, 2)
        + 5*Power(w, 2)*Power(x, 2) + 4*v*Power(x, 3))
    + Power(w, 3)*(Power(w, 6) - 3*v*Power(w, 4)*x
        - Power(w, 2)*Power(x, 3) + v*(Power(v, 2) - 2*x)*Power(x, 3))
    - 3*u*Power(w, 2)*(v*Power(w, 6) - Power(w, 4)*Power(x, 2)
        + Power(v, 4)*Power(x, 3) + Power(x, 5)
        - 3*Power(v, 2)*x*(Power(w, 4) + Power(x, 3)))
    + Power(u, 4)*w*x*(-12*Power(v, 3)*Power(w, 2)*x
        + 3*Power(w, 4)*x + 3*Power(v, 4)*Power(x, 2)
        - 23*v*Power(w, 2)*Power(x, 2) - 4*Power(x, 4)
        + Power(v, 2)*(3*Power(w, 4) + 4*Power(x, 3)))
    + Power(u, 2)*w*(3*Power(w, 6)*x + 3*Power(v, 5)*Power(x, 3)
        - 4*Power(w, 2)*Power(x, 4) + 9*v*Power(x, 2)*(-2*Power(w, 4)
          + Power(x, 3)) + 3*Power(v, 2)*(Power(w, 6)
          + 2*Power(w, 2)*Power(x, 3)) - 3*Power(v, 3)*(3*Power(w, 4)*x
          + 4*Power(x, 4)))
    + Power(u, 3)*(-(Power(v, 6)*Power(x, 3))
        + 7*Power(w, 4)*Power(x, 3) + Power(x, 6)
        - Power(v, 3)*(Power(w, 6) + 8*Power(w, 2)*Power(x, 3))
        + Power(v, 4)*(3*Power(w, 4)*x + 5*Power(x, 4))
        + v*(-6*Power(w, 6)*x + 4*Power(w, 2)*Power(x, 4))
        + 3*Power(v, 2)*(9*Power(w, 4)*Power(x, 2) - 2*Power(x, 5)));
  b[0][0] = b[0][0]/den;

  b[0][1] = Power(u, 7)*Power(x, 3)*(2*Power(w, 2) - 3*v*x)
    - Power(u, 6)*w*Power(x, 2)*(6*v*Power(w, 2) - 10*Power(v, 2)*x
        + Power(x, 2))
    - Power(w, 3)*(-(Power(v, 3)*Power(w, 2)*x)
        + 2*Power(w, 4)*x - Power(v, 4)*Power(x, 2) + Power(x, 4)
        + Power(v, 2)*(Power(w, 4) + Power(x, 3)))
    + Power(u, 5)*x*(-6*Power(v, 3)*Power(w, 2)*x + 6*Power(w, 4)*x
        - 5*Power(v, 4)*Power(x, 2) - 10*v*Power(w, 2)*Power(x, 2)
        + 2*Power(v, 2)*(3*Power(w, 4) + Power(x, 3)))
    - Power(u, 2)*w*(-3*Power(v, 5)*Power(w, 2)*x
        + 3*Power(v, 2)*Power(w, 4)*x - 3*Power(v, 6)*Power(x, 2)
        - 12*Power(v, 3)*Power(w, 2)*Power(x, 2)
        + 3*Power(x, 2)*(Power(w, 4) - Power(x, 3))
        + 3*Power(v, 4)*(Power(w, 4) + 3*Power(x, 3))
        + v*(6*Power(w, 6) + 5*Power(w, 2)*Power(x, 3)))
    - Power(u, 4)*w*(Power(v, 4)*Power(w, 2)*x - 6*Power(v, 5)*Power(x, 2)
        - 12*Power(v, 2)*Power(w, 2)*Power(x, 2) + Power(w, 2)*Power(x, 3)
        + 2*Power(v, 3)*(Power(w, 4) - Power(x, 3))
        + 4*v*(3*Power(w, 4)*x - Power(x, 4)))
    + u*Power(w, 2)*(2*Power(w, 6) - 3*Power(v, 4)*Power(w, 2)*x
        - 3*Power(v, 5)*Power(x, 2) - 3*Power(v, 2)*Power(w, 2)*Power(x, 2)
        + 3*Power(v, 3)*(Power(w, 4) + 2*Power(x, 3))
        + v*(4*Power(w, 4)*x + 3*Power(x, 4)))
    - Power(u, 3)*(Power(v, 6)*Power(w, 2)*x - 6*Power(w, 6)*x
        + Power(v, 7)*Power(x, 2) + 15*Power(v, 4)*Power(w, 2)*Power(x, 2)
        + 5*Power(w, 2)*Power(x, 4)
        - Power(v, 5)*(Power(w, 4) + 4*Power(x, 3))
        - 2*Power(v, 2)*(3*Power(w, 6) + 4*Power(w, 2)*Power(x, 3))
        + Power(v, 3)*(-2*Power(w, 4)*x + 2*Power(x, 4))
        + v*(3*Power(w, 4)*Power(x, 2) + 2*Power(x, 5)));
  b[1][0] = b[0][1] = b[0][1]/den;

  b[0][2] = Power(u, 8)*w*Power(x, 3)
    - Power(u, 7)*Power(x, 2)*(3*v*Power(w, 2) - 4*Power(v, 2)*x
        + Power(x, 2))
    + Power(u, 6)*w*x*(3*Power(v, 2)*Power(w, 2) - 4*Power(v, 3)*x
        + 3*Power(w, 2)*x - 7*v*Power(x, 2))
    - Power(w, 3)*(2*v*Power(w, 4) - 2*Power(v, 2)*Power(w, 2)*x
        - 2*Power(v, 3)*Power(x, 2) + Power(w, 2)*Power(x, 2)
        + 3*v*Power(x, 3))
    + u*Power(w, 2)*(-6*Power(v, 3)*Power(w, 2)*x
        - 6*Power(v, 4)*Power(x, 2) + v*Power(w, 2)*Power(x, 2)
        - 3*Power(x, 4) + 3*Power(v, 2)*(2*Power(w, 4) + 5*Power(x, 3)))
    + Power(u, 4)*w*(-9*Power(v, 3)*Power(w, 2)*x
        + 5*Power(v, 4)*Power(x, 2) - 12*v*Power(w, 2)*Power(x, 2)
        + 3*x*(Power(w, 4) - 2*Power(x, 3))
        + Power(v, 2)*(3*Power(w, 4) + 16*Power(x, 3)))
    + Power(u, 3)*(-2*Power(v, 5)*Power(w, 2)*x
        - 2*Power(v, 6)*Power(x, 2) - 14*Power(v, 3)*Power(w, 2)*Power(x, 2)
        + Power(x, 2)*(Power(w, 4) + Power(x, 3))
        + Power(v, 4)*(2*Power(w, 4) + 9*Power(x, 3))
        - v*(3*Power(w, 6) + 2*Power(w, 2)*Power(x, 3))
        + 3*Power(v, 2)*(5*Power(w, 4)*x - 3*Power(x, 4)))
    + Power(u, 5)*(Power(v, 4)*Power(w, 2)*x + Power(v, 5)*Power(x, 2)
        + 15*Power(v, 2)*Power(w, 2)*Power(x, 2) + Power(w, 2)*Power(x, 3)
        - Power(v, 3)*(Power(w, 4) + 12*Power(x, 3))
        + v*(-6*Power(w, 4)*x + 8*Power(x, 4)))
    + Power(u, 2)*w*(Power(w, 6) + 6*Power(v, 4)*Power(w, 2)*x
        + 6*Power(v, 5)*Power(x, 2)
        + 8*Power(v, 2)*Power(w, 2)*Power(x, 2)
        - 2*Power(w, 2)*Power(x, 3)
        - 3*Power(v, 3)*(2*Power(w, 4) + 7*Power(x, 3))
        + v*(-7*Power(w, 4)*x + 12*Power(x, 4)));
  b[2][0] = b[0][2] = b[0][2]/den;

  b[0][3] = -(Power(u, 6)*w*Power(x, 3))
    + Power(u, 5)*Power(x, 2)*(3*v*Power(w, 2) - 4*Power(v, 2)*x
        + Power(x, 2))
    + Power(u, 4)*w*x*(-3*Power(v, 2)*Power(w, 2)
        + 4*Power(v, 3)*x - 3*Power(w, 2)*x + 5*v*Power(x, 2))
    - Power(w, 3)*(Power(w, 4) - v*Power(w, 2)*x
        + Power(x, 2)*(-Power(v, 2) + x))
    - u*Power(w, 2)*(3*Power(v, 2)*Power(w, 2)*x
        + 3*Power(v, 3)*Power(x, 2) + Power(w, 2)*Power(x, 2)
        - 3*v*(Power(w, 4) + 2*Power(x, 3)))
    + 3*Power(u, 2)*w*(Power(v, 3)*Power(w, 2)*x - Power(w, 4)*x
        + Power(v, 4)*Power(x, 2) + 2*v*Power(w, 2)*Power(x, 2)
        + Power(x, 4) - Power(v, 2)*(Power(w, 4) + 3*Power(x, 3)))
    - Power(u, 3)*(Power(v, 4)*Power(w, 2)*x + Power(v, 5)*Power(x, 2)
        + 9*Power(v, 2)*Power(w, 2)*Power(x, 2) + Power(w, 2)*Power(x, 3)
        - Power(v, 3)*(Power(w, 4) + 4*Power(x, 3))
        + v*(-6*Power(w, 4)*x + 3*Power(x, 4)));
  b[3][0] = b[0][3] = b[0][3]/den;

  b[1][1] = 3*Power(u, 8)*w*Power(x, 3)
    - 3*Power(u, 7)*Power(x, 2)*(4*v*Power(w, 2) + Power(v, 2)*x
        + Power(x, 2))
    + Power(u, 6)*w*x*(12*Power(v, 2)*Power(w, 2)
        + 16*Power(v, 3)*x + 12*Power(w, 2)*x + 15*v*Power(x, 2))
    + Power(w, 3)*(Power(v, 4)*Power(w, 2) + Power(v, 5)*x
        + 4*Power(v, 2)*Power(w, 2)*x + 2*Power(w, 2)*Power(x, 2)
        - v*Power(x, 3))
    - u*Power(w, 2)*(3*Power(v, 5)*Power(w, 2)
        + 3*Power(v, 6)*x + 16*Power(v, 3)*Power(w, 2)*x + 8*Power(w, 4)*x
        - 3*Power(v, 4)*Power(x, 2) + 6*v*Power(w, 2)*Power(x, 2)
        - 3*Power(x, 4) + Power(v, 2)*(4*Power(w, 4) - 6*Power(x, 3)))
    + Power(u, 2)*w*(3*Power(v, 6)*Power(w, 2) + 4*Power(w, 6)
        + 3*Power(v, 7)*x + 27*Power(v, 4)*Power(w, 2)*x
        + 28*v*Power(w, 4)*x - 6*Power(v, 5)*Power(x, 2)
        + 8*Power(v, 2)*Power(w, 2)*Power(x, 2) - 3*Power(w, 2)*Power(x, 3)
        + 6*Power(v, 3)*(2*Power(w, 4) - Power(x, 3)))
    - Power(u, 5)*(16*Power(v, 4)*Power(w, 2)*x + 24*v*Power(w, 4)*x
        + 5*Power(v, 5)*Power(x, 2) + 36*Power(v, 2)*Power(w, 2)*Power(x, 2)
        + 17*Power(w, 2)*Power(x, 3)
        + 4*Power(v, 3)*(Power(w, 4) + Power(x, 3)))
    + Power(u, 4)*w*(4*Power(v, 5)*Power(w, 2) + 7*Power(v, 6)*x
        + 44*Power(v, 3)*Power(w, 2)*x + 9*Power(v, 4)*Power(x, 2)
        + 36*v*Power(w, 2)*Power(x, 2)
        + 4*Power(v, 2)*(3*Power(w, 4) + 2*Power(x, 3))
        + 4*(3*Power(w, 4)*x + Power(x, 4)))
    - Power(u, 3)*(Power(v, 7)*Power(w, 2) + 12*v*Power(w, 6)
        + Power(v, 8)*x + 22*Power(v, 5)*Power(w, 2)*x
        - 3*Power(v, 6)*Power(x, 2) + 8*Power(v, 3)*Power(w, 2)*Power(x, 2)
        + 16*Power(w, 4)*Power(x, 2) + Power(x, 5)
        + Power(v, 4)*(12*Power(w, 4) - Power(x, 3))
        + 2*Power(v, 2)*(24*Power(w, 4)*x + Power(x, 4)));
  b[1][1] =  b[1][1]/den;

  b[1][2] = -6*Power(u, 8)*v*w*Power(x, 2) + Power(u, 9)*Power(x, 3)
    + 2*Power(u, 7)*x*(3*Power(v, 2)*Power(w, 2) + 2*Power(v, 3)*x
        + 3*Power(w, 2)*x)
    + Power(w, 3)*(2*Power(v, 3)*Power(w, 2)
        + 2*Power(v, 4)*x + 4*v*Power(w, 2)*x - Power(v, 2)*Power(x, 2)
        - Power(x, 3))
    - Power(u, 6)*w*(2*Power(v, 3)*Power(w, 2)
        + 5*Power(v, 4)*x + 12*v*Power(w, 2)*x - 6*Power(v, 2)*Power(x, 2)
        + 5*Power(x, 3))
    + u*Power(w, 2)*(-6*Power(v, 4)*Power(w, 2)
        - 4*v*Power(w, 4) - 6*Power(v, 5)*x - 16*Power(v, 2)*Power(w, 2)*x
        + 9*Power(v, 3)*Power(x, 2) + 2*Power(w, 2)*Power(x, 2)
        + 6*v*Power(x, 3))
    + Power(u, 4)*w*(Power(v, 4)*Power(w, 2)
        - 6*v*Power(w, 4) + 7*Power(v, 5)*x + 9*Power(v, 2)*Power(w, 2)*x
        + 16*Power(v, 3)*Power(x, 2) - Power(w, 2)*Power(x, 2)
        + 10*v*Power(x, 3))
    + Power(u, 2)*w*(6*Power(v, 5)*Power(w, 2)
        + 6*Power(v, 6)*x + 29*Power(v, 3)*Power(w, 2)*x - 2*Power(w, 4)*x
        - 15*Power(v, 4)*Power(x, 2) - 4*v*Power(w, 2)*Power(x, 2)
        + 3*Power(x, 4) + Power(v, 2)*(11*Power(w, 4) - 3*Power(x, 3)))
    + Power(u, 5)*(Power(v, 5)*Power(w, 2) + Power(v, 6)*x
        + 2*Power(v, 3)*Power(w, 2)*x - 11*Power(v, 4)*Power(x, 2)
        - 9*v*Power(w, 2)*Power(x, 2)
        + Power(v, 2)*(6*Power(w, 4) - 2*Power(x, 3))
        + 2*x*(3*Power(w, 4) + Power(x, 3)))
    - Power(u, 3)*(2*Power(v, 6)*Power(w, 2) - 2*Power(w, 6)
        + 2*Power(v, 7)*x + 25*Power(v, 4)*Power(w, 2)*x
        - 7*Power(v, 5)*Power(x, 2) + 3*Power(v, 2)*Power(w, 2)*Power(x, 2)
        + 5*Power(w, 2)*Power(x, 3)
        + Power(v, 3)*(9*Power(w, 4) + 2*Power(x, 3))
        + v*(4*Power(w, 4)*x + 3*Power(x, 4)));
  b[2][1] = b[1][2] = b[1][2]/den;

  b[1][3] = Power(v, 2)*Power(w, 5) + Power(v, 3)*Power(w, 3)*x
    + 2*Power(w, 5)*x + 6*Power(u, 6)*v*w*Power(x, 2)
    - Power(u, 7)*Power(x, 3)
    - 2*Power(u, 5)*x*(3*Power(v, 2)*Power(w, 2)
        + 2*Power(v, 3)*x + 3*Power(w, 2)*x + v*Power(x, 2))
    - u*Power(w, 2)*(3*Power(v, 3)*Power(w, 2) + 2*Power(w, 4)
        + 3*Power(v, 4)*x + 8*v*Power(w, 2)*x - 3*Power(v, 2)*Power(x, 2)
        - 3*Power(x, 3))
    + Power(u, 4)*w*(2*Power(v, 3)*Power(w, 2)
        + 5*Power(v, 4)*x + 12*v*Power(w, 2)*x + 6*Power(v, 2)*Power(x, 2)
        + 4*Power(x, 3))
    + Power(u, 2)*w*(3*Power(v, 4)*Power(w, 2)
        + 6*v*Power(w, 4) + 3*Power(v, 5)*x + 15*Power(v, 2)*Power(w, 2)*x
        - 6*Power(v, 3)*Power(x, 2) + Power(w, 2)*Power(x, 2)
        - 3*v*Power(x, 3))
    - Power(u, 3)*(Power(v, 5)*Power(w, 2) + 6*Power(v, 2)*Power(w, 4)
        + Power(v, 6)*x + 14*Power(v, 3)*Power(w, 2)*x + 6*Power(w, 4)*x
        - 3*Power(v, 4)*Power(x, 2) + 3*v*Power(w, 2)*Power(x, 2)
        + Power(x, 4));
  b[3][1] = b[1][3] = b[1][3]/den;

  b[2][2] = -3*Power(u, 9)*v*Power(x, 2) + 3*Power(u, 8)*w*x*(Power(v, 2) + x)
    + 4*v*Power(w, 3)*(v*Power(w, 2) + Power(v, 2)*x - Power(x, 2))
    + 3*Power(u, 6)*w*(Power(v, 2)*Power(w, 2) - 3*Power(v, 3)*x
        + Power(w, 2)*x - 6*v*Power(x, 2))
    - 3*u*Power(w, 2)*(4*Power(v, 3)*Power(w, 2) + 4*Power(v, 4)*x
        - 8*Power(v, 2)*Power(x, 2) + Power(x, 3))
    - Power(u, 7)*(Power(v, 3)*Power(w, 2) + Power(v, 4)*x
        + 6*v*Power(w, 2)*x - 15*Power(v, 2)*Power(x, 2) + 3*Power(x, 3))
    + Power(u, 2)*w*(12*Power(v, 4)*Power(w, 2) - 4*v*Power(w, 4)
        + 12*Power(v, 5)*x + 8*Power(v, 2)*Power(w, 2)*x
        - 36*Power(v, 3)*Power(x, 2) + 2*Power(w, 2)*Power(x, 2)
        + 15*v*Power(x, 3))
    + Power(u, 4)*(-12*Power(v, 3)*Power(w, 3)
        + Power(w, 5) - 11*v*Power(w, 3)*x + 42*Power(v, 2)*w*Power(x, 2)
        - 7*w*Power(x, 3))
    + Power(u, 3)*(-4*Power(v, 5)*Power(w, 2)
        - 4*Power(v, 6)*x - 12*Power(v, 3)*Power(w, 2)*x
        + 16*Power(v, 4)*Power(x, 2) - 18*v*Power(w, 2)*Power(x, 2)
        + Power(x, 4) + Power(v, 2)*(12*Power(w, 4) - 13*Power(x, 3)))
    + Power(u, 5)*(4*Power(v, 4)*Power(w, 2) + 4*Power(v, 5)*x
        + 21*Power(v, 2)*Power(w, 2)*x - 26*Power(v, 3)*Power(x, 2)
        + 3*Power(w, 2)*Power(x, 2) - 3*v*(Power(w, 4) - 4*Power(x, 3)));
  b[2][2] = b[2][2]/den;

  b[2][3] = 3*Power(u, 7)*v*Power(x, 2) - 3*Power(u, 6)*w*x*(Power(v, 2) + x)
    - 3*u*v*Power(w, 2)*(2*v*Power(w, 2) + 2*Power(v, 2)*x - 3*Power(x, 2))
    + Power(w, 3)*(2*v*Power(w, 2) + 2*Power(v, 2)*x - Power(x, 2))
    - 3*Power(u, 4)*w*(Power(v, 2)*Power(w, 2) - Power(v, 3)*x
        + Power(w, 2)*x - 4*v*Power(x, 2))
    + Power(u, 5)*(Power(v, 3)*Power(w, 2) + Power(v, 4)*x
        + 6*v*Power(w, 2)*x - 9*Power(v, 2)*Power(x, 2) + 2*Power(x, 3))
    + Power(u, 2)*w*(6*Power(v, 3)*Power(w, 2) - Power(w, 4)
        + 6*Power(v, 4)*x + 5*v*Power(w, 2)*x - 15*Power(v, 2)*Power(x, 2)
        + 3*Power(x, 3))
    - Power(u, 3)*(2*Power(v, 4)*Power(w, 2) - 3*v*Power(w, 4)
        + 2*Power(v, 5)*x + 9*Power(v, 2)*Power(w, 2)*x
        - 7*Power(v, 3)*Power(x, 2) + 3*Power(w, 2)*Power(x, 2)
        + 4*v*Power(x, 3));
  b[3][2] = b[2][3] = b[2][3]/den;

  b[3][3] = Power(w, 5) + v*Power(w, 3)*x - 3*Power(u, 5)*v*Power(x, 2)
    + 3*Power(u, 4)*w*x*(Power(v, 2) + x)
    - 3*u*Power(w, 2)*(v*Power(w, 2) + Power(v, 2)*x - Power(x, 2))
    + 3*Power(u, 2)*w*(Power(v, 2)*Power(w, 2) + Power(v, 3)*x
        + Power(w, 2)*x - 2*v*Power(x, 2))
    - Power(u, 3)*(Power(v, 3)*Power(w, 2) + Power(v, 4)*x
        + 6*v*Power(w, 2)*x - 3*Power(v, 2)*Power(x, 2) + Power(x, 3));
  b[3][3] = b[3][3]/den;

#ifdef TIMING
  TOC(4, time_compute_fhb)
#endif
}
// -----------------------------------------------------------------
