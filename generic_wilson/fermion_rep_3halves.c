/*
 Create link field 'link' in rep of fermions from field 'linkf'
 in fundamental rep
*/

/* fermion_rep_3halves.c:  construct Isospin=3/2 irrep

			 -- -- --
			|  |  |  |
			 -- -- --

   Basis states are:

	[000 001 011 111]
*/

#include "generic_wilson_includes.h"

#if (NCOL != 2)
	#error "Wrong version of fermion_rep!"
#endif
#if (DIMF != 4)
	#error "Wrong version of fermion_rep!"
#endif
#if (FREP != 3/2)
	#error "Wrong version of fermion_rep!"
#endif

/*  "hard coded" generation of I=3/2 group elements  */
void make_fermion_rep_matrix(su3_matrix_f *a, su3_matrix *b){

  complex u00,u01,u10,u11,u00sq,u01sq,u10sq,u11sq,t,s;
  double root3;

  root3=sqrt(3.);

  /* copy fund_rep matrix and reunitarize "on the fly" */
  u00=a->e[0][0];
  u01=a->e[0][1];
  CONJG(u01,t);  CNEGATE(t,u10);
  CONJG(u00,u11);

  /* squares */
  CMUL(u00,u00,u00sq);
  CMUL(u01,u01,u01sq);
  CMUL(u10,u10,u10sq);
  CMUL(u11,u11,u11sq);

  /* now the I=3/2 matrix */
  CMUL(u00,u00sq,b->e[0][0]);      /* u00^3                        */
  CMUL(u01,u00sq,t);
  CMULREAL(t,root3,b->e[0][1]);    /* sqrt(3)*u00^2*u01            */
  CMUL(u00,u01sq,t);
  CMULREAL(t,root3,b->e[0][2]);    /* sqrt(3)*u01^2*u00            */
  CMUL(u01,u01sq,b->e[0][3]);      /* u01^3                        */
  CMUL(u10,u00sq,t);
  CMULREAL(t,root3,b->e[1][0]);    /* sqrt(3)*u00^2*u10            */
  CMUL(u10,u00,t);
  CMUL(t,u01,s);
  CMULREAL(s,2.,t);
  CMUL(u11,u00sq,s);
  CADD(s,t,b->e[1][1]);            /* u00^2*u11 + 2*u10*u00*u01    */
  CMUL(u11,u00,t);
  CMUL(t,u01,s);
  CMULREAL(s,2.,t);
  CMUL(u10,u01sq,s);
  CADD(s,t,b->e[1][2]);            /* u01^2*u10 + 2*u11*u00*u01    */
  CMUL(u11,u01sq,t);
  CMULREAL(t,root3,b->e[1][3]);    /* sqrt(3)*u01^2*u11            */
  CMUL(u00,u10sq,t);
  CMULREAL(t,root3,b->e[2][0]);    /* sqrt(3)*u10^2*u00            */
  CMUL(u11,u00,t);
  CMUL(t,u10,s);
  CMULREAL(s,2.,t);
  CMUL(u01,u10sq,s);
  CADD(s,t,b->e[2][1]);            /* u10^2*u01 + 2*u11*u00*u10    */
  CMUL(u10,u11,t);
  CMUL(t,u01,s);
  CMULREAL(s,2.,t);
  CMUL(u00,u11sq,s);
  CADD(s,t,b->e[2][2]);            /* u11^2*u00 + 2*u10*u11*u01    */
  CMUL(u01,u11sq,t);
  CMULREAL(t,root3,b->e[2][3]);    /* sqrt(3)*u11^2*u01            */
  CMUL(u10,u10sq,b->e[3][0]);      /* u10^3                        */
  CMUL(u11,u10sq,t);
  CMULREAL(t,root3,b->e[3][1]);    /* sqrt(3)*u10^2*u11            */
  CMUL(u10,u11sq,t);
  CMULREAL(t,root3,b->e[3][2]);    /* sqrt(3)*u11^2*u10            */
  CMUL(u11,u11sq,b->e[3][3]);      /* u11^3                        */

} /* END make_fermion_rep_matrix()	*/



/*** make_fermion_rep_driv *******************************************
YS Sept 07

Generates parametric derivatives of boundary links using
the Leibniz rule, for use in measuring 1/g^2(L).

Minimilly adapted from make_fermion_rep_matrix.
The code assumes that the boundary matrices are diagonal

*/

#ifdef SF

#define MULT_BY_THREE_FACTORS(a,f_one,f_two) {   \
  complex temp;                                  \
  temp.imag = (a).real*( 2*(f_one) + (f_two) );  \
  temp.real = -(a).imag*( 2*(f_one) + (f_two) ); \
  (a) = temp;                                    \
  }

void make_fermion_rep_driv(su3_matrix_f *a, su3_matrix *b,
                           Real *factor, int sign) {
  int i,j;
  complex u00,u11,u00sq,u11sq;
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
  CMUL(u00,u00,u00sq);
  CMUL(u11,u11,u11sq);

  /* now the I=3/2 parametric derivative */
  CMUL(u00sq,u00,b->e[0][0]);
  MULT_BY_THREE_FACTORS(b->e[0][0], sfactor[0], sfactor[0]);
  CMUL(u00sq,u11,b->e[1][1]);
  MULT_BY_THREE_FACTORS(b->e[1][1], sfactor[0], sfactor[1]);
  CMUL(u11sq,u00,b->e[2][2]);
  MULT_BY_THREE_FACTORS(b->e[2][2], sfactor[1], sfactor[0]);
  CMUL(u11sq,u11,b->e[3][3]);
  MULT_BY_THREE_FACTORS(b->e[3][3], sfactor[1], sfactor[1]);

}	/* END make_fermion_rep_driv */

#endif /* SF */

