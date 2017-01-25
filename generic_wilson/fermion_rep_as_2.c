/*
 Create link field 'link' in rep of fermions from field 'linkf'
 in fundamental rep
*/

/* fermion_rep_as_2.c :  construct two-index antisymmetric rep

			 --
			|  |
			 --
			|  |
			 --

   Basis states are:

	[01 02 03 ... 12 13 ... 23 ...]
*/

#include "generic_wilson_includes.h"

#if (DIMF != NCOL*(NCOL-1)/2)
	#error "Wrong version of fermion_rep!"
#endif
#if (FREP != antisymmetric2)
        #error "Wrong version of fermion_rep!"
#endif

void make_fermion_rep_matrix(matrix_f *a, matrix *b){

int i,j,ij,k,l,kl; /* ij and kl are the compound indices of the
				fermion rep  */
complex x,y;

ij=0;
for(i=0;i<NCOL;i++) for(j=i+1;j<NCOL;j++){
	kl=0;
	for(k=0;k<NCOL;k++) for(l=k+1;l<NCOL;l++){
		CMUL(a->e[i][k], a->e[j][l], x);
		CMUL(a->e[i][l], a->e[j][k], y);
		CSUB(x, y, b->e[ij][kl]);
		kl++;
	}
	ij++;
}
}	/* END make_fermion_rep_matrix() */


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

void make_fermion_rep_driv(matrix_f *a, matrix *b,
                           Real *factor, int sign) {

int i,j,ij,k,l,kl; /* ij and kl are the compound indices of the
				fermion rep  */
complex x,y;
Real sfactor[NCOL];

for(i=0;i<NCOL;i++) {
  sfactor[i] = (Real)sign*factor[i];
}

ij=0;
for(i=0;i<NCOL;i++) for(j=i+1;j<NCOL;j++){
	kl=0;
	for(k=0;k<NCOL;k++) for(l=k+1;l<NCOL;l++){
		CMUL(a->e[i][k], a->e[j][l], x);
		CMUL(a->e[i][l], a->e[j][k], y);
		CSUB(x, y, b->e[ij][kl]);
                MULT_BY_FACTORS(b->e[ij][kl], sfactor[i], sfactor[j]);
		kl++;
	}
	ij++;
}
}	/* END make_fermion_rep_driv */

#endif /* SF */
