/*
 Create link field 'link' in rep of fermions from field 'linkf'
 in fundamental rep
*/

/* fermion_rep_sym_2.c :  construct two-index symmetric rep

			 -- --
			|  |  |
			 -- --

   Basis states are:

	[00 01 02 03 ... 11 12 13 ... 22 23 ... 33 ...]
*/

#include "generic_wilson_includes.h"

#if (DIMF != NCOL*(NCOL+1)/2)
	#error "Wrong version of fermion_rep!"
#endif
#if (FREP != symmetric2)
        #error "Wrong version of fermion_rep!"
#endif

void make_fermion_rep_matrix(matrix_f *a, matrix *b){

int i,j,ij,k,l,kl; /* ij and kl are the compound indices of the
				fermion rep  */
complex x,y;
double root2;

root2=sqrt(2.);

ij=0;
for(i=0;i<NCOL;i++) {
   kl=0;
   for(k=0;k<NCOL;k++){					/* j=i */
      CMUL(a->e[i][k], a->e[i][k], b->e[ij][kl]);       /* l=k */
      kl++;
      for(l=k+1;l<NCOL;l++){			        /* l>k */
	 CMUL(a->e[i][k], a->e[i][l], x);
	 CMULREAL(x, root2, b->e[ij][kl]);
	 kl++;
      }
   }
   ij++;
   for(j=i+1;j<NCOL;j++){				/* j>i */
      kl=0;
      for(k=0;k<NCOL;k++){
         CMUL(a->e[i][k], a->e[j][k], x);	        /* l=k */
         CMULREAL(x, root2, b->e[ij][kl]);
         kl++;
         for(l=k+1;l<NCOL;l++){		        /* l>k */
            CMUL(a->e[i][k], a->e[j][l], x);
            CMUL(a->e[i][l], a->e[j][k], y);
            CADD(x, y, b->e[ij][kl]);
            kl++;
         }
      }
      ij++;
   }
}
}
