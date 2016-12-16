/*
  Create link field 'link' in rep of fermions from field 'linkf'
  in fundamental rep


  fermion_rep_2halves.c:  HARD CODED construction of isospin=1

			 -- --
			|  |  |
			 -- --

   Basis states are:

	[00 01 11]
*/

#include "generic_wilson_includes.h"

#if (NCOL != 2)
	#error "Wrong version of fermion_rep!"
#endif
#if (DIMF != 3)
	#error "Wrong version of fermion_rep!"
#endif
#if (FREP != symmetric2)
	#error "Wrong version of fermion_rep!"
#endif

void make_fermion_rep_matrix(su3_matrix_f *a, su3_matrix *b){

  complex u00,u01,u10,u11,t;
  double root2=sqrt(2.);

  /* copy fund_rep matrix and reunitarize "on the fly" */
  u00=a->e[0][0];
  u01=a->e[0][1];
  CONJG(u01,t);  CNEGATE(t,u10);
  CONJG(u00,u11);

  /* now the I=1 matrix */
  CMUL(u00,u00,b->e[0][0]);        /* u00^2                        */
  CMUL(u00,u01,t);
  CMULREAL(t,root2,b->e[0][1]);    /* sqrt(2)*u00*u01              */
  CMUL(u01,u01,b->e[0][2]);        /* u01^2                        */
  CMUL(u00,u10,t);
  CMULREAL(t,root2,b->e[1][0]);    /* sqrt(2)*u00*u10              */
  CMUL(u00,u11,t);
  CMUL(u10,u01,b->e[1][1]);
  CSUM(b->e[1][1],t);              /* u00*u11 + u10*u01            */
  CMUL(u01,u11,t);
  CMULREAL(t,root2,b->e[1][2]);    /* sqrt(2)*u01*u11              */
  CMUL(u10,u10,b->e[2][0]);        /* u10^2                        */
  CMUL(u10,u11,t);
  CMULREAL(t,root2,b->e[2][1]);    /* sqrt(2)*u10*u11              */
  CMUL(u11,u11,b->e[2][2]);        /* u11^2                        */


} /* END make_fermion_rep_matrix()	*/


/*** make_fermion_rep_driv *******************************************
YS Sept 07

Generates parametric derivatives of boundary links using
the Leibniz rule, for use in measuring 1/g^2(L).

Minimilly adapted from make_fermion_rep_matrix.
The code assumes that the boundary matrices are diagonal

*/

#ifdef SF

#define MULT_BY_FACTORS(a,f_one,f_two) {       \
  complex temp;                                \
  temp.imag = (a).real*( (f_one) + (f_two) );  \
  temp.real = -(a).imag*( (f_one) + (f_two) ); \
  (a) = temp;                                  \
  }

void make_fermion_rep_driv(su3_matrix_f *a, su3_matrix *b,
                           Real *factor, int sign) {
  int i,j;
  complex u00,u11;
  Real sfactor[NCOL];

  for(i=0;i<NCOL;i++) {
    sfactor[i] = (Real)sign*factor[i];
  }

  for(i=0;i<DIMF;i++)for(j=0;j<DIMF;j++){
    b->e[i][j]=cmplx(0.0,0.0);
  }

  /* copy diagonal elements of fund_rep matrix */
  u00=a->e[0][0];
  CONJG(u00,u11);

  /* now the I=1 parametric derivative */
  CMUL(u00,u00,b->e[0][0]);
  MULT_BY_FACTORS(b->e[0][0], sfactor[0], sfactor[0]);
  CMUL(u00,u11,b->e[1][1]);
  MULT_BY_FACTORS(b->e[1][1], sfactor[0], sfactor[1]);
  CMUL(u11,u11,b->e[2][2]);
  MULT_BY_FACTORS(b->e[2][2], sfactor[1], sfactor[1]);


} /* END make_fermion_rep_driv */

#endif /* SF */
