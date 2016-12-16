/***************** nhyp.c - SU(2) version ******************************/

#if (NCOL != 2)
        #error "Wrong NCOL in nhyp.c!"
#endif

/* compute_fhb
  This routine computes the Cayley Hamilton coefficients for the inverse sqrt
  such that 1/sqrt(Q) = f[0] Id for SU(2) (i.e., f[1] = 0 always)
  The b[][] matrix contains the derivatives of f wrt to the traces of Q
  b[i][j] = d f[i] / d c[j] with c[j] = 1/(j+1) trace Q^(j+1)
  The flag compute_b switches on/off the computation of these derivatives
*/

#ifndef NHYP_DEBUG
void compute_fhb( su3_matrix_f *Q, Real *f, Real b[NCOL][NCOL], int compute_b )
#else
void compute_fhb( su3_matrix_f *Omega, su3_matrix_f *Q,
                  Real *f, Real b[NCOL][NCOL], int compute_b )
#endif

{
  int i,j;
  Real c0;
  /* this is easy. We already have Q= Omega^\dagger Omega so */

  c0=Q->e[0][0].real+Q->e[1][1].real;
  f[0]=sqrt(2.0/c0);
  f[1]=0.;

  if(compute_b){
    for(i=0;i<2;i++) for(j=0;j<2;j++)	b[i][j]=0.0;
    b[0][0]= - sqrt(0.5/c0/c0/c0);
  }
}
