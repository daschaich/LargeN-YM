/*
 Translate force term, calculated in higher rep, to fund rep
*/

/* fund_rep_force_3halves.c :  force given in Isospin=3/2 irrep
*
*    		         -- -- --
*			|  |  |  |
*			 -- -- --
*
*   Basis states are:
*
*	[000 001 011 111]
*
*/

#include "cl_dyn_includes.h"

#if (NCOL != 2)
	#error "Wrong version of fermion_rep!"
#endif
#if (DIMF != 4)
	#error "Wrong version of fermion_rep!"
#endif
#if (FREP != 3/2)
        #error "Wrong version of fermion_rep!"
#endif

/*  "hard coded" generation of force in fundamental irrep

        / 3/2    0    0    0 \          /   0 sqrt(3)   0       0 \
        |   0  1/2    0    0 |          |   0       0   2       0 |
  T_3 = |   0    0 -1/2    0 |    T_+ = |   0       0   0 sqrt(3) |
        \   0    0    0 -3/2 /          \   0       0   0       0 /

  pi-dot = 1/2 sum_{a=1}^3 sigma_a tr[T_a(A-A^dag)]

  where T_a, A are in the I=3/2 irrep.

*/

void fund_rep_force(su3_matrix_f *force, su3_matrix *a){

  Real diag;
  complex s,t,u;
  double root3=sqrt(3.);

  /* off-diagonal generators                               */

  CADD(a->e[0][1],a->e[2][3],t);
  CONJG(t,s);
  CADD(a->e[1][0],a->e[3][2],t);
  CSUB(t,s,u);
  CMULREAL(u,root3,u);
  CONJG(a->e[1][2],t);
  CSUB(a->e[2][1],t,s);
  CMULREAL(s,2.,t);
  CADD(u,t,s);
/* now fix overall normalization of off-diagonal terms */
  CMULREAL(s,0.5,force->e[1][0]);

/*
now do e[0][1] = -e[1][0]^* .
needed since make_anti_hermitian uses both triangular parts.
*/
	force->e[0][1].real = -force->e[1][0].real;
        force->e[0][1].imag = force->e[1][0].imag;

  /* diagonal generator (sigma_3) */

  diag = 1.5*(a->e[0][0].imag - a->e[3][3].imag)
        +0.5*(a->e[1][1].imag - a->e[2][2].imag);
  force->e[0][0].imag = diag;
  force->e[1][1].imag = -diag;
/* set diag.real to zero - this doesn't affect result,
only to have force look consistent.
*/
  force->e[0][0].real=0;
  force->e[1][1].real=0;

} /* END fund_rep_force()	*/
