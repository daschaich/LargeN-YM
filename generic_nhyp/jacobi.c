// -----------------------------------------------------------------
#include <math.h>
#include "../include/complex.h"
#include "../include/jacobi.h"

//#define DEBUG
#ifdef DEBUG
#include <stdio.h>
#endif

#define sign(a) ( a<0 ? -1.0 : 1.0)
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Allocate space for an NxN Matrix
Matrix AllocateMatrix(int N) {
  register int i;
  Matrix A;

  A.N = N;
  A.M = malloc(N * sizeof(double_complex *));
  for (i = 0; i < N; i++)
    A.M[i] = malloc(N * sizeof(double_complex));
  return A;
}

// Deallocate the space of Matrix A
void deAllocate(Matrix *A) {
  register int i, N = A->N;
  for (i = 0; i < N; i++)
    free(A->M[i]);
  free(A->M);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Return the Frobenius norm of the matrix A,
//   sqrt(Sum_ij |A_ij|^2)
Real FrobeniusNorm(Matrix *A) {
  register int i, j, N = A->N;
  register Real tmp = 0.0;
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      tmp += A->M[i][j].real * A->M[i][j].real
           + A->M[i][j].imag * A->M[i][j].imag;
    }
  }
  return sqrt(tmp);
}

// Return only sqrt(Sum_ij(i!=j) |A_ij|^2)
Real OffDiag(Matrix *A) {
  register int i, j, N = A->N;
  register Real tmp = 0.0;
  for (i = 0; i < N; i++) {
    for (j = i + 1; j < N; j++) {
      tmp += A->M[i][j].real * A->M[i][j].real
           + A->M[i][j].imag * A->M[i][j].imag;
      tmp += A->M[j][i].real * A->M[j][i].real
           + A->M[j][i].imag * A->M[j][i].imag;
    }
  }
  return sqrt(tmp);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Set B = Adag
void HermitianConj(Matrix *A, Matrix *B) {
  register int i, j, N = A->N;
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      B->M[i][j].real =  A->M[j][i].real;
      B->M[i][j].imag = -A->M[j][i].imag;
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute the Schur decomposition of the hermitian matrix A,
// setting the (p, q) element to zero, with p < q
// Reference: Page 446 of "Matrix Computations" by Golub & Van Loan, 1989
void HermitianSchur(Matrix *A, int p, int q, h_schur *schur) {
  register Real b, t, off, c, s, Apq_r, Apq_i;

  schur->p = p;
  schur->q = q;
  Apq_r = A->M[p][q].real;
  Apq_i = A->M[p][q].imag;
  off = sqrt(Apq_r * Apq_r + Apq_i * Apq_i);
  if (off > DEPS_SCHUR) {
    off = 1.0 / off;
    b = 0.5 * (A->M[q][q].real - A->M[p][p].real) * off;
    t = sign(b) / (fabs(b) + sqrt(b * b + 1));
    c = 1.0 / sqrt(1 + t * t);
    schur->c = c;
    s = c * t * off;
    schur->s.real = Apq_r * s;
    schur->s.imag = Apq_i * s;
  }
  else {
    schur->c = 1.0;
    schur->s.real = 0.0;
    schur->s.imag = 0.0;
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Multiply from the right with the Schur rotation
void rightMultiplySchur(Matrix *A, h_schur *schur) {
  register int i, N = A->N, p = schur->p, q = schur->q;
  register Real c = schur->c, s_r, s_i, Tp_r, Tp_i, Tq_r, Tq_i;

  s_r  = schur->s.real;
  s_i  = schur->s.imag;

  for (i = 0; i < N; i++) {
    Tp_r = A->M[i][p].real;
    Tp_i = A->M[i][p].imag;
    Tq_r = A->M[i][q].real;
    Tq_i = A->M[i][q].imag;

    A->M[i][p].real = Tp_r * c  - Tq_r * s_r - Tq_i * s_i;
    A->M[i][p].imag = Tp_i * c  + Tq_r * s_i - Tq_i * s_r;

    A->M[i][q].real = Tq_r * c  + Tp_r * s_r - Tp_i * s_i;
    A->M[i][q].imag = Tq_i * c  + Tp_r * s_i + Tp_i * s_r;
  }
}

// Multiply from the right with the (hermitian conjugated) Schur rotation
void leftMultiplySchur(Matrix *A, h_schur *schur) {
  register int i, N = A->N, p = schur->p, q = schur->q;
  register Real c = schur->c, s_r, s_i, Tp_r, Tp_i, Tq_r, Tq_i;

  s_r  = schur->s.real;
  s_i  = schur->s.imag;

  for (i = 0; i < N; i++) {
    Tp_r = A->M[p][i].real;
    Tp_i = A->M[p][i].imag;
    Tq_r = A->M[q][i].real;
    Tq_i = A->M[q][i].imag;

    A->M[p][i].real = Tp_r * c  - Tq_r * s_r + Tq_i * s_i;
    A->M[p][i].imag = Tp_i * c  - Tq_r * s_i - Tq_i * s_r;

    A->M[q][i].real = Tq_r * c  + Tp_r * s_r + Tp_i * s_i;
    A->M[q][i].imag = Tq_i * c  - Tp_r * s_i + Tp_i * s_r;
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute the Matrix V that diagonalizes Vdag.A.V
void Jacobi(Matrix *A, Matrix *V, Real tolerance) {
  register int i, j, N = A->N, iter = 0;
  register Real eps, off, old_off;
  h_schur S;

  // Initialize V to the unit matrix
  for (i=0;i<N;i++) {
    V->M[i][i].real = 1.0;
    V->M[i][i].imag = 0.0;
    for (j = i + 1; j < N; j++) {
      V->M[i][j].real = 0.0;
      V->M[i][j].imag = 0.0;
      V->M[j][i].real = 0.0;
      V->M[j][i].imag = 0.0;
    }
  }

  eps = FrobeniusNorm(A) * tolerance;
#ifdef DEBUG
  printf("In jacobi::Jacobi -- convergence eps: %g\n",eps);
  printf("%i off(A): %g\n", iter, OffDiag(A));
#endif
  old_off  = 1.0e+32;
  off = OffDiag(A);
  while(off > eps
     && fabs(old_off - off) > eps
     && iter < MAX_JACOBI_ITERS) {
    iter++;
    for (i = 0; i < N; i++) {
      for (j = i + 1; j < N; j++) {
        HermitianSchur(A, i, j, &S);
        rightMultiplySchur(A, &S);
        leftMultiplySchur(A, &S);
        rightMultiplySchur(V, &S);
      }
    }
    old_off = off;
    off = OffDiag(A);
#ifdef DEBUG
    printf("%i off(A): %g\n", iter, off);
#endif
  }
  if (iter < JACOBI_HIST_MAX)
    jacobi_hist[iter]++;
  if (iter == MAX_JACOBI_ITERS)
    node0_printf("Jacobi didn't converge in %d steps\n", MAX_JACOBI_ITERS);
}
// -----------------------------------------------------------------
