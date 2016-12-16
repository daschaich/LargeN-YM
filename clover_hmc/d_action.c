/*************** d_action.c ****************************************/
/* MIMD version 6 */

/* Measure total action, as needed by the hybrid Monte Carlo algorithm.
   When this routine is called the conjugate gradient should already
   have been run on the even sites, so that the vector psi contains
   (M_adjoint*M)^(-1) * chi.
*/

#include "cl_dyn_includes.h"
Real ahmat_mag_sq(anti_hermitmat *pt);

double d_action(){
double d_hmom_action(),d_fermion_action();
double ssplaq,g_action,h_action,f_action;
#ifndef IMP
double stplaq;
#ifdef BETA_FREP
double ssplaq_frep,stplaq_frep,g_action_frep;
#endif /* BETA_FREP */
#endif

#ifdef IMP
    gauge_action(&ssplaq);
    g_action = beta*ssplaq/(Real)NCOL;
#else
/* d_plaquette returns average ss and st plaqs       */
    d_plaquette(&ssplaq, &stplaq);
    ssplaq=1-ssplaq/(Real)NCOL;
    stplaq=1-stplaq/(Real)NCOL;
#ifndef SF
    g_action = beta*3*nx*ny*nz*nt*(ssplaq+stplaq);
#else
    g_action = beta*3*nx*ny*nz*((nt-1)*ssplaq+nt*stplaq);
#endif /* not SF  */
#ifdef BETA_FREP
    d_plaquette_frep(&ssplaq_frep, &stplaq_frep);
    ssplaq_frep=1-ssplaq_frep/(Real)DIMF;
    stplaq_frep=1-stplaq_frep/(Real)DIMF;
#ifndef SF
    g_action_frep = beta_frep*3*nx*ny*nz*nt*(ssplaq_frep+stplaq_frep);
#else
    g_action_frep = beta_frep*3*nx*ny*nz*((nt-1)*ssplaq_frep+nt*stplaq_frep);
#endif /* not SF  */
#endif /* BETA_FREP */
#endif /* IMP     */

    h_action = d_hmom_action();
    f_action = d_fermion_action();
if(this_node==0)printf("D_ACTION: g,h,f = %e  %e  %e  %e\n",
g_action,h_action,f_action,
(g_action+h_action+f_action));
/*  return(g_action+h_action); */
#ifndef BETA_FREP
    return(g_action+h_action+f_action);
#else
    return(g_action+g_action_frep+h_action+f_action);
#endif
}

/* fermion contribution to the action */
double d_fermion_action() {
register int i;
register site *s;
 double sum;
 sum= (double)0.0;
#ifndef LU
   FORALLSITESDOMAIN(i,s){
#else
     FOREVENSITESDOMAIN(i,s){
#endif
	if(num_masses==1)
      	 sum += (double)wvec_rdot( &(s->psi[0]), &(s->chi[0]) );
	else{
       /* recall that level 0 is the slow invert and 1 is the fast (shifted) one */
	sum += (double)wvec_rdot( &(s->chi[0]), &(s->chi[0]) );
       	sum+= shift*shift*(double)wvec_rdot( &(s->psi[0]), &(s->chi[0]) );
	sum += (double)wvec_rdot( &(s->psi[1]), &(s->chi[1]) );
	}
     }
    g_doublesum(&sum);
/**{Real xxx ; xxx = sum; g_floatsum( &xxx ); sum = xxx;}**/
    return(sum);
}

/* gauge momentum contribution to the action */
double d_hmom_action() {
register int i,dir;
register site *s;
double sum;

    sum= (double)0.0;
    for(dir=XUP;dir<=TUP;dir++){
        FORALLDYNLINKS(i,s,dir){
            sum += (double)ahmat_mag_sq( &(s->mom[dir]) );
        }
    }
    g_doublesum( &sum );
/**{Real xxx ; xxx = sum; g_floatsum( &xxx ); sum = xxx;}**/
    return(sum);
}

/* magnitude squared of an antihermition matrix */
Real ahmat_mag_sq(anti_hermitmat *pt){
register Real x,sum;
    x = pt->m00im; sum  = 0.5*x*x;
    x = pt->m11im; sum += 0.5*x*x;
    x = pt->m01.real; sum += x*x;
    x = pt->m01.imag; sum += x*x;
#if (NCOL>2)
    x = pt->m22im; sum += 0.5*x*x;
    x = pt->m02.real; sum += x*x;
    x = pt->m02.imag; sum += x*x;
    x = pt->m12.real; sum += x*x;
    x = pt->m12.imag; sum += x*x;
#if (NCOL>3)
    x = pt->m33im; sum += 0.5*x*x;
    x = pt->m03.real; sum += x*x;
    x = pt->m03.imag; sum += x*x;
    x = pt->m13.real; sum += x*x;
    x = pt->m13.imag; sum += x*x;
    x = pt->m23.real; sum += x*x;
    x = pt->m23.imag; sum += x*x;
#endif
#endif
    return(sum);
}

/* copy a gauge field - an array of four su3_matrices */
/* SF: don't care that we copy also t=0 boundary links */
void gauge_field_copy_f(field_offset src,field_offset dest) {
register int i,dir,src2,dest2;
register site *s;
    FORALLSITES(i,s){
	src2=src; dest2=dest;
        for(dir=XUP;dir<=TUP; dir++){
	    su3mat_copy_f( (su3_matrix_f *)F_PT(s,src2),
		(su3_matrix_f *)F_PT(s,dest2) );
	    src2 += sizeof(su3_matrix_f);
	    dest2 += sizeof(su3_matrix_f);
	}
    }
}
void gauge_field_copy(field_offset src,field_offset dest) {
register int i,dir,src2,dest2;
register site *s;
    FORALLSITES(i,s){
	src2=src; dest2=dest;
        for(dir=XUP;dir<=TUP; dir++){
	    su3mat_copy( (su3_matrix *)F_PT(s,src2),
		(su3_matrix *)F_PT(s,dest2) );
	    src2 += sizeof(su3_matrix);
	    dest2 += sizeof(su3_matrix);
	}
    }
}
