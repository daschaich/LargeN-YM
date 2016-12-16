/************************** path.c *******************************/
/* MIMD version 6 */
/* An arbitrary path walker */
/*  uses tempmat1.e[0][0] to accumulate, and ordinary gathers*/

/* Added _f suffixes throughout for fundamental rep operations. 
	Work exclusively with linkf field, not link.  -bqs */

#include "cl_dyn_includes.h"

void path(int *dir,int *sign,int length)
{
register int i;
register site *s;
msg_tag *mtag0;
int j;


/* j=0 */
	if(sign[0]>0)  {
	    mtag0 = start_gather_site( F_OFFSET(linkf[dir[0]]), sizeof(su3_matrix_f),
		OPP_DIR(dir[0]), EVENANDODD, gen_pt[0] );
	    wait_gather(mtag0);

	      FORALLSITES(i,s){
	      su3mat_copy_f((su3_matrix_f *)(gen_pt[0][i]),(su3_matrix_f *)&(s->tempmat1) );
	      }

	    cleanup_gather(mtag0);
	}

	if(sign[0]<0) { 
	      FORALLSITES(i,s){
	      su3_adjoint_f(&(s->linkf[dir[0]]),(su3_matrix_f *)&(s->tempmat1) );
	      }
	}


	for(j=1;j<length;j++) {
	if(sign[j] > 0) {

	      FORALLSITES(i,s){
		mult_su3_nn_f( (su3_matrix_f *)&(s->tempmat1),&(s->linkf[dir[j]]),
		    (su3_matrix_f *)&(s->tempmat2) );
	      }
	    mtag0 = start_gather_site( F_OFFSET(tempmat2), sizeof(su3_matrix_f),
		OPP_DIR(dir[j]), EVENANDODD, gen_pt[0] );
	    wait_gather(mtag0);

	      FORALLSITES(i,s){
	      su3mat_copy_f((su3_matrix_f *)(gen_pt[0][i]),(su3_matrix_f *)&(s->tempmat1) ); 
	      }
	    cleanup_gather(mtag0);
	}

	if(sign[j] < 0) {
	    mtag0 = start_gather_site( F_OFFSET(tempmat1), sizeof(su3_matrix_f),
		dir[j], EVENANDODD, gen_pt[1] );
	    wait_gather(mtag0);

	    FORALLSITES(i,s){
		mult_su3_na_f((su3_matrix_f *)(gen_pt[1][i]),
			&(s->linkf[dir[j]]), (su3_matrix_f *)&(s->tempmat2) );
	    }

	    FORALLSITES(i,s){
	      su3mat_copy_f((su3_matrix_f *)&(s->tempmat2),(su3_matrix_f *)&(s->tempmat1) );
	    }
	    cleanup_gather(mtag0);
	}

      }



} /* path */

