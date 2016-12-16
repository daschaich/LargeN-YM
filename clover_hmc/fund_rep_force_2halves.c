/*
 Translate force term, calculated in higher rep, to fund rep
*/

/* fund_rep_force_2halves.c :   HARD CODED construction of isospin=1

			 -- --
			|  |  |
			 -- --

   Basis states are:

	[00 01 11]
*/

#include "cl_dyn_includes.h"

#if (NCOL != 2)
	#error "Wrong version of fermion_rep!"
#endif
#if (DIMF != 3)
	#error "Wrong version of fermion_rep!"
#endif
#if (FREP != symmetric2)
	#error "Wrong version of fermion_rep!"
#endif


/*  "hard coded" generation of force in fundamental irrep

        / +1   0   0 \          / 0  sqrt(2)       0 \
  T_3 = |  0   0   0 |    T_+ = | 0        0 sqrt(2) |
        \  0   0  -1 /          \ 0        0       0 /

  pi-dot = 1/2 sum_{a=1}^3 sigma_a tr[T_a(B-B^dag)]

  where T_a, B are in the I=1 irrep.

*/

void fund_rep_force(su3_matrix_f *force, su3_matrix *a){

  Real diag;
  complex s,t;
  double root2=sqrt(2.);

  /* off-diagonal generators                               */

  CADD(a->e[0][1],a->e[1][2],t);
  CONJG(t,s);
  CADD(a->e[1][0],a->e[2][1],t);
  CSUB(t,s,force->e[1][0]);
  CMULREAL(force->e[1][0],root2/2.,force->e[1][0]);

/*
now do e[0][1] = -e[1][0]^* .
needed since make_anti_hermitian uses both triangular parts.
*/
	force->e[0][1].real = -force->e[1][0].real;
        force->e[0][1].imag = force->e[1][0].imag;

  /* diagonal generator (sigma_3) */

  diag = a->e[0][0].imag - a->e[2][2].imag;
  force->e[0][0].imag = diag;
  force->e[1][1].imag = -diag;

/* set diag.real to zero - this doesn't affect result,
only to have force look consistent.
resulting force is already traceless and anti-hermitian.
*/
  force->e[0][0].real=0;
  force->e[1][1].real=0;

} /* END fund_rep_force()	*/
