// -----------------------------------------------------------------
// Reunitarization for arbitrary NCOL
// 1) Do a singular value decomposition (SVD) using LAPACK for arbitrary NCOL
//      in --> L.S.Rdag
// 2) Reconstruct out = L.Rdag, setting vector S=(1, 1, ..., 1)
// Both L and R are unitary NCOLxNCOL matrices
#include "generic_includes.h"

#define TOLERANCE 0.0001
#define MAXERRCOUNT 100
//#define UNIDEBUG

Real max_deviation;
double av_deviation;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int check_deviation(Real deviation) {
  if (max_deviation < deviation)
    max_deviation = deviation;
  av_deviation += deviation * deviation;

  if (deviation < TOLERANCE)
    return 0;
  else
    return 1;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Reunitarize using LAPACK SVD
int reunit(matrix *c) {
  register int i;
  char A = 'A';     // Ask LAPACK for all singular values
  int row, col, Npt = NCOL, stat = 0, Nwork = 3 * NCOL, err = 0;
  Real dev;
  matrix lmat, rdagmat;

  // Convert c to column-major double array used by LAPACK
  for (row = 0; row < NCOL; row++) {
    for (col = 0; col < NCOL; col++) {
      i = 2 * (col * NCOL + row);
      store[i] = c->e[row][col].real;
      store[i + 1] = c->e[row][col].imag;
    }
  }

  // Compute singular value decomposition of c
  zgesvd_(&A, &A, &Npt, &Npt, store, &Npt, junk, left, &Npt, right, &Npt,
          work, &Nwork, Rwork, &stat);

  // Move the results back into matrix structures
  clear_mat(&lmat);
  clear_mat(&rdagmat);
  for (row = 0; row < NCOL; row++) {
    for (col = 0; col < NCOL; col++) {
      i = 2 * (col * NCOL + row);
      lmat.e[row][col].real = left[i];
      lmat.e[row][col].imag = left[i + 1];
      rdagmat.e[row][col].real = right[i];
      rdagmat.e[row][col].imag = right[i + 1];
    }
  }

  // Now u = l.rdag (throwing out singular values in junk)
  mult_nn(&lmat, &rdagmat, c);

  // Check the unitarity of the result
  dev = check_unit(c);
  err = check_deviation(dev);
  if (err)
    dump_mat(c);
  return err;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void reunitarize() {
  register int i, dir;
  register site *s;
  register matrix *mat;
  int errcount = 0, errors;

  max_deviation = 0.0;
  av_deviation = 0.0;

  FORALLSITES(i, s) {
    FORALLUPDIR(dir) {
      mat = &(s->link[dir]);
      errors = reunit(mat);
      errcount += errors;
      if (errors) {
        printf("Unitarity problem above (node %d, site %d, tolerance %.4g)\n",
               mynode(), i, TOLERANCE);
      }
      if (errcount > MAXERRCOUNT) {
        printf("Unitarity error count exceeded: %d\n", errcount);
        fflush(stdout);
        terminate(1);
      }
    }
  }

  av_deviation = sqrt(av_deviation / (double)nt);
#ifdef UNIDEBUG
  printf("Deviation from unitarity on node %d: max %.4g, ave %.4g\n",
         mynode(), max_deviation, av_deviation);
#endif
  if (max_deviation > TOLERANCE) {
    printf("reunitarize: Node %d unitarity problem, maximum deviation %.4g\n",
           mynode(), max_deviation);
    errcount++;
    if (errcount > MAXERRCOUNT) {
      printf("Unitarity error count exceeded\n");
      terminate(1);
    }
  }
}
// -----------------------------------------------------------------
