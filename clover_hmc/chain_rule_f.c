/********** chain_rule_f2.c ****************************/
/* MIMD version 7 */
/* YS 8/08	  */

/* chain_rule_f.c :  dummy, force already given in fund rep
   chain_rule_f.c is just as trivial as fund_rep_force_f.c,
   but we need to be consistent with chain_rule(..) prototype.

			 --
			|  |
			 --
*/

#include "cl_dyn_includes.h"

#if (DIMF != NCOL)
	#error "Wrong version of fund_rep_force!"
#endif
#if (FREP != fundamental)
	#error "Wrong version of fund_rep_force!"
#endif

void chain_rule(su3_matrix_f *sigmaf, su3_matrix *sigma,
                su3_matrix_f *gaugelinkf){

  int i,j;

  for(i=0;i<NCOL;i++) {
    for(j=0;j<NCOL;j++){
      sigmaf->e[i][j]=sigma->e[i][j];
    } /* end loop over j */
  } /* end loop over i */
} /* END chain-rule() */


