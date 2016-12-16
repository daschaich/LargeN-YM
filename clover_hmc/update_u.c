/*************** update_u.c ************************************/
/* MIMD version 6 */
/* Same as wilson_dynamical/update_u.c */

/* update the link matrices					*
*  								*
*  Go to sixth order in the exponential of the momentum		*
*  matrices, since unitarity seems important.			*
*  Evaluation is done as:					*
*	exp(H) * U = ( 1 + H + H^2/2 + H^3/3 ...)*U		*
*	= U + H*( U + (H/2)*( U + (H/3)*( ... )))		*
*								*
*/

#include "cl_dyn_includes.h"

void update_u(Real eps){

register int i,dir;
register site *s;
su3_matrix_f *link,temp1,temp2,htemp;
register Real t2,t3,t4,t5,t6;

/* Temporary by-hand optimization until pgcc compiler bug is fixed */
t2 = eps/2.0;
t3 = eps/3.0;
t4 = eps/4.0;
t5 = eps/5.0;
t6 = eps/6.0;

    for(dir=XUP; dir <=TUP; dir++){
        FORALLDYNLINKS(i,s,dir){
	    uncompress_anti_hermitian( &(s->mom[dir]) , &htemp );
	    link = &(s->linkf[dir]);
	    mult_su3_nn_f(&htemp,link,&temp1);
            /**scalar_mult_add_su3_matrix(link,&temp1,eps/6.0,&temp2);**/
scalar_mult_add_su3_matrix_f(link,&temp1,t6,&temp2);
	    mult_su3_nn_f(&htemp,&temp2,&temp1);
            /**scalar_mult_add_su3_matrix(link,&temp1,eps/5.0,&temp2);**/
scalar_mult_add_su3_matrix_f(link,&temp1,t5,&temp2);
	    mult_su3_nn_f(&htemp,&temp2,&temp1);
            /**scalar_mult_add_su3_matrix(link,&temp1,eps/4.0,&temp2);**/
scalar_mult_add_su3_matrix_f(link,&temp1,t4,&temp2);
	    mult_su3_nn_f(&htemp,&temp2,&temp1);
	    /**scalar_mult_add_su3_matrix(link,&temp1,eps/3.0,&temp2);**/
scalar_mult_add_su3_matrix_f(link,&temp1,t3,&temp2);
	    mult_su3_nn_f(&htemp,&temp2,&temp1);
	    /**scalar_mult_add_su3_matrix(link,&temp1,eps/2.0,&temp2);**/
scalar_mult_add_su3_matrix_f(link,&temp1,t2,&temp2);
	    mult_su3_nn_f(&htemp,&temp2,&temp1);
	    scalar_mult_add_su3_matrix_f(link,&temp1,eps    ,&temp2);
	    su3mat_copy_f(&temp2,link);
        }
    }

/**dtime += dclock();
if(this_node==0)printf("LINK_UPDATE: time = %e  mflops = %e\n",
dtime, (double)(5616.0*volume/(1.0e6*dtime*numnodes())) );**/
} /* update_u */
