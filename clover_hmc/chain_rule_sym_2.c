/*
 Translate force term, calculated in higher rep, to fund rep,
 via the chain rule.  Closely follows make_fermion_rep_matrix
*/

/* fund_rep_force_sym_2.c :  force given in two-index symmetric rep

			 -- --
			|  |  |
			 -- --

   Basis states are:

	[00 01 02 03 ... 11 12 13 ... 22 23 ...]
*/

#include "cl_dyn_includes.h"

#if (DIMF != NCOL*(NCOL+1)/2)
	#error "Wrong version of fund_rep_force!"
#endif
#if (FREP != symmetric2)
        #error "Wrong version of fund_rep_force!"
#endif

void chain_rule(su3_matrix_f *sigmaf, su3_matrix *sigma,
                su3_matrix_f *gaugelinkf){

  int i,j,ij,k,l,kl; /* ij and kl are the compound indices of the
				fermion rep  */
  complex x,y;
  double root2=sqrt(2.);

  clear_su3mat_f(sigmaf);

  ij=0;
  for(i=0;i<NCOL;i++) {
    kl=0;                       		/* j=i */
    for(k=0;k<NCOL;k++){
      CMUL(gaugelinkf->e[i][k], sigma->e[kl][ij], x);        /* l=k */
      CMULREAL(x, 2., y);
      CSUM(sigmaf->e[k][i], y);
      kl++;
      for(l=k+1;l<NCOL;l++){			        /* l>k */
        CMUL(gaugelinkf->e[i][k], sigma->e[kl][ij], x);
        CMULREAL(x, root2, y);
        CSUM(sigmaf->e[l][i], y);
        CMUL(gaugelinkf->e[i][l], sigma->e[kl][ij], x);
        CMULREAL(x, root2, y);
        CSUM(sigmaf->e[k][i], y);
        kl++;
      } /* end loop over l */
    } /* end loop over k */
    ij++;
    for(j=i+1;j<NCOL;j++){			 /* j>i */
      kl=0;
      for(k=0;k<NCOL;k++){
        CMUL(gaugelinkf->e[i][k], sigma->e[kl][ij], x);      /* l=k */
        CMULREAL(x, root2, y);
        CSUM(sigmaf->e[k][j], y);
        CMUL(gaugelinkf->e[j][k], sigma->e[kl][ij], x);
        CMULREAL(x, root2, y);
        CSUM(sigmaf->e[k][i], y);
        kl++;
        for(l=k+1;l<NCOL;l++){		                /* l>k */
          CMUL(gaugelinkf->e[i][k], sigma->e[kl][ij], y);
          CSUM(sigmaf->e[l][j], y);
          CMUL(gaugelinkf->e[j][k], sigma->e[kl][ij], y);
          CSUM(sigmaf->e[l][i], y);
          CMUL(gaugelinkf->e[i][l], sigma->e[kl][ij], y);
          CSUM(sigmaf->e[k][j], y);
          CMUL(gaugelinkf->e[j][l], sigma->e[kl][ij], y);
          CSUM(sigmaf->e[k][i], y);
          kl++;
        } /* end loop over l */
      } /* end loop over k */
      ij++;
    } /* end loop over j */
  } /* end loop over i */
} /* END chain-rule() */


