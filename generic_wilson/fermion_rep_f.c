/*
 Create link field 'link' in rep of fermions from field 'linkf'
 in fundamental rep
*/

/* fermion_rep_f.c :  dummy, transfers fundamental rep from linkf to link

			 --
			|  |
			 --

   Basis states are:

	[0 1 2 3 4 ... ]
*/

#include "generic_wilson_includes.h"

#if (DIMF != NCOL)
	#error "Wrong version of fermion_rep!"
#endif
#if (FREP != fundamental)
        #error "Wrong version of fermion_rep!"
#endif


void make_fermion_rep_matrix(su3_matrix_f *a, su3_matrix *b){

int i,j;

for(i=0;i<NCOL;i++) for(j=0;j<NCOL;j++)
	b->e[i][j] = a->e[i][j];

}	/* END make_fermion_rep_matrix() */


/*** make_fermion_rep_driv *******************************************
YS Sept 07

Generates parametric derivatives of boundary links using
the Leibniz rule, for use in measuring 1/g^2(L).

Minimilly adapted from make_fermion_rep_matrix.
The code assumes that the boundary matrices are diagonal

*/

#ifdef SF

void make_fermion_rep_driv(su3_matrix_f *a, su3_matrix *b,
                           Real *factor, int sign) {
  int i,j;

  for(i=0;i<NCOL;i++)for(j=0;j<NCOL;j++){
    b->e[i][j]=cmplx(0.0,0.0);
  }

  for(i=0;i<NCOL;i++) {
    b->e[i][i].imag = (Real)sign * factor[i] * a->e[i][i].real;
    b->e[i][i].real = -(Real)sign * factor[i] * a->e[i][i].imag;
  }
}	/* END make_fermion_rep_driv */

#endif /* SF */
