/*
 Translate force term, calculated in higher rep, to fund rep,
 via the chain rule.  Closely follows make_fermion_rep_matrix
*/

/* fund_rep_force_as_2.c :  force given in two-index antisymmetric rep

                         --
                        |  |
                         --
                        |  |
                         --

   Basis states are:

        [01 02 03 ... 12 13 ... 23 ...]

*/

#include "cl_dyn_includes.h"

#if (DIMF != NCOL*(NCOL-1)/2)
	#error "Wrong version of fund_rep_force!"
#endif
#if (FREP != antisymmetric2)
        #error "Wrong version of fermion_rep!"
#endif

void chain_rule(matrix_f *sigmaf, matrix *sigma,
                matrix_f *gaugelinkf){

  int i,j,ij,k,l,kl; /* ij and kl are the compound indices of the
				fermion rep  */
  complex y;

  clear_mat_f(sigmaf);

  ij=0;
  for(i=0;i<NCOL;i++) {
    for(j=i+1;j<NCOL;j++){			 /* j>i */
      kl=0;
      for(k=0;k<NCOL;k++){
        for(l=k+1;l<NCOL;l++){		                /* l>k */
          CMUL(gaugelinkf->e[i][k], sigma->e[kl][ij], y);
          CSUM(sigmaf->e[l][j], y);
          CMUL(gaugelinkf->e[j][k], sigma->e[kl][ij], y);
          CSUB(sigmaf->e[l][i], y, sigmaf->e[l][i]);
          CMUL(gaugelinkf->e[i][l], sigma->e[kl][ij], y);
          CSUB(sigmaf->e[k][j], y, sigmaf->e[k][j]);
          CMUL(gaugelinkf->e[j][l], sigma->e[kl][ij], y);
          CSUM(sigmaf->e[k][i], y);
          kl++;
        } /* end loop over l */
      } /* end loop over k */
      ij++;
    } /* end loop over j */
  } /* end loop over i */
} /* END chain-rule() */


