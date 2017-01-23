// -----------------------------------------------------------------
#ifndef _JACOBI_H
#define _JACOBI_H
#define DEPS_SCHUR 1.110223e-16
// MAX_JACOBI_ITERS is now defined in defines.h
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// The Schur rotation of a hermitian matrix
typedef struct {
  double c;
  double_complex s;
  int p;
  int q;
} h_schur;

typedef struct {
  int N;
  double_complex **M;
} Matrix;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
Matrix AllocateMatrix(int N);
void deAllocate(Matrix *A);

Real FrobeniusNorm(Matrix *A);
Real OffDiag(Matrix *A);

void HermitianConj(Matrix *A, Matrix *B);
void HermitianSchur(Matrix *A, int p, int q, h_schur *schur);
void rightMultiplySchur(Matrix *A, h_schur *schur);
void leftMultiplySchur(Matrix *A, h_schur *schur);
void Jacobi(Matrix *A, Matrix *V, Real tolerance);

#endif
// -----------------------------------------------------------------
