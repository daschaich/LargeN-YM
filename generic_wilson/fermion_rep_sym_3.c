/********** fermion_rep_sym_3.c ****************************/

/*  NOT FINISHED!  */

/* MIMD version 7 */
/* BS 1/07	*/

/*
 Create link field 'link' in rep of fermions from field 'linkf'
 in fundamental rep
*/

/* fermion_rep_sym_3.c :  construct three-index symmetric rep

			 -- -- --
			|  |  |  |
			 -- -- --

   Definition of vector in this rep in terms of fundamental indices:

	v_a = [v_111 v_112 v_113 v_114 ... v_122 v_123 v_124 ... v_133 v_134 ...]
*/

#include "generic_wilson_includes.h"

#if (DIMF != NCOL*(NCOL+1)*(NCOL+2)/6)
	#error "Wrong version of fermion_rep!"
#endif

void make_fermion_rep_matrix(su3_matrix_f *a, su3_matrix *b){

int i,j,r,ijr,k,l,s,kls; /* ij and kl are the compound indices of the
				fermion rep  */
complex x;

ijr=0;
for(i=0;i<NCOL;i++) for(j=i;j<NCOL;j++) for(r=j;r<NCOL;r++){
	kls=0;
	for(k=0;k<NCOL;k++){
		CMUL(a->e[i][k], a->e[j][k], b->e[ij][kl]);
		kl++;
		for(l=k+1;l<NCOL;l++){
			CMUL(a->e[i][k], a->e[j][l], b->e[ij][kl]);
			CMUL(a->e[i][l], a->e[j][k], x);
			CSUM(b->e[ij][kl], x);
			kl++;
		}
	}
	ij++;
}
