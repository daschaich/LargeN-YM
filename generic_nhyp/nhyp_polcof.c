/***************** nhyp.c - SU(4) version ******************************/
#include "jacobi.c"

#include "Power.h"
#include "../generic/nrutil.h"

/* compute_fhb
  This routine computes the Cayley Hamilton coefficients for the inverse sqrt
  such that 1/sqrt(Q) = f[0] Id + f[1] Q + f[2] Q^2 + f[3] Q^3
  using a general Vandermonde formula from Numerical Recipes.

It doesn't compute the b's! If compute_b=1, the code dies.

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
  
  /*
  double sg0,sg1,sg2;
  double u1, p, den;
  double u12,u02;
  double u0;
*/
  double g[NCOL],q[NCOL],ff[NCOL];
  
  
  
  /*  Matrix Qj, Vj; moved to lattice.h   */
  
#ifdef TIMING
  TIC(4)
#endif
    
    /**** Find eigenvalues of Q *******************************************/
    
    /*  Qj = AllocateMatrix(NCOL);
	Vj = AllocateMatrix(NCOL);  moved to setup.c */
    
    
    for(i=0;i<NCOL;i++)
      for(j=0;j<NCOL;j++){
	Qj.M[i][j].real=Q->e[i][j].real;
	Qj.M[i][j].imag=Q->e[i][j].imag;
      }
  /*  
      for(i=0;i<NCOL;i++){
      for(j=0;j<NCOL;j++)printf("%e %e   ",Qj.M[i][j].real,Qj.M[i][j].imag);
      printf("\n");
      }
  */
  
  
  Jacobi(&Qj,&Vj,TOL_JACOBI);
  for(i=0;i<NCOL;i++)g[i]=Qj.M[i][i].real;

  /*
  for(i=0;i<NCOL;i++)printf("GG %e\n ",g[i]);
  */
 
#ifdef NHYP_DEBUG
  for(i=0;i<NCOL;i++){
    if(Qj.M[i][i].real < 0.0){
      printf("NHYP_DEBUG_FATAL: negative eigenvalue\n  = %.12e\n",
	     Qj.M[i][i].real);
      DUMP_STUFF
	terminate(1);
    }
  }
#endif
 

  for(i=0;i<NCOL;i++) q[i]=1.0/sqrt(g[i]);
 
 
  polcof(g,q,NCOL-1,ff); /* Note N.R. convention, vector is 0,1,2,...N, NOT c-like! */

  /*
  for(i=0;i<NCOL;i++){printf("POLF %d %e\n",i,ff[i]);}
  for(i=0;i<NCOL;i++){
    u0= ff[0];
    for(j=1;j<NCOL;j++) {
      u0 += ff[j]*Power(g[i],j);
    }
    printf("CHECKV %i %e %e %e \n",i,u0,q[i],u0-q[i]);
  }
  printf("END NR\n");
  */
  for(i=0;i<NCOL;i++)f[i]=ff[i]; /*ff in  NR is in double...*/
  
  /* compare to SU3 analytic case 
     
  sg0 = sqrt(g[0]);
  sg1 = sqrt(g[1]);
  sg2 = sqrt(g[2]);
  
  u0 = sg0+sg1+sg2;
  u1 = sg0*sg1+sg0*sg2+sg1*sg2;
  p  = sg0*sg1*sg2;
  
  u02=u0*u0;
  u12=u1*u1;
  
  den=p*(u0*u1-p);
  f[0]=(-p*(u02 + u1) + u0*u12)/den;
  f[1]=(-p + u0* (2*u1-u02))/den;
  f[2]=u0/den;
  
  for(i=0;i<NCOL;i++){printf("%d %e %e\n",i,f[i],g[i] );}
  for(i=0;i<NCOL;i++){
  u0= ff[0];
  for(j=1;j<NCOL;j++) {
  u0 += f[j]*Power(g[i],j);
  }
  printf("CHECKV %i %e %e %e \n",i,u0,q[i],u0-q[i]);
  }
  printf("END ANALYTIC\n");
  
  end su3 analytic */
  
  /* most of the time, the coefficients are all we need */
  if (compute_b!=0){
    node0_printf("Can't compute b's\n");
    exit(1);
  }
  else{
    
#ifdef TIMING
    TOC(4,time_compute_fhb)
#endif
      }
  return;
}
