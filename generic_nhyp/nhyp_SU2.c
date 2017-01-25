// -----------------------------------------------------------------
// Helper routines for SU(2) nHYP blocking

#if NCOL != 2
  #error "Wrong NCOL in nhyp.c!"
#endif

// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute the Cayley Hamilton coefficients for the inverse square root
//   1 / sqrt(Q) = f[0] I (although f[1] exists, it is zero)
// The b[][] matrix contains the derivatives of f wrt to the traces of Q
//   b[i][j] = d f[i] / d c[j] with c[j] = 1/(j+1) trace Q^(j+1)
// The flag compute_b switches on/off the computation of these derivatives
#ifndef NHYP_DEBUG
void compute_fhb( matrix_f *Q, Real *f, Real b[NCOL][NCOL], int compute_b )
#else
void compute_fhb( matrix_f *Omega, matrix_f *Q,
                  Real *f, Real b[NCOL][NCOL], int compute_b )
#endif

{
  int i,j;
  Real c0;
  /* this is easy. We already have Q= Omega^\dagger Omega so */

  c0 = Q->e[0][0].real + Q->e[1][1].real;
  f[0] = sqrt(2.0 / c0);
  f[1] = 0.;

  if(compute_b) {
    for (i=0;i<2;i++) for (j=0;j<2;j++) b[i][j]=0.0;
    b[0][0]= - sqrt(0.5/c0/c0/c0);
  }
}
// -----------------------------------------------------------------
