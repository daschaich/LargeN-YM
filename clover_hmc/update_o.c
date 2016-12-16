/*
 Update lattice.
Omelian integrator with Hasenbusch


Use Urbach, Jansen,Schindler, Wenger multiscale, CPC 174 (2006) 87


Note that for the final accept/reject, we already have a good solution to the CG; the last
update was of the momenta.

The is in a .h file:
nsteps[num_masses+1]
We  read in the number of levels in the input file.
Usually, ``levels''=number of masses+1, where the number of masses is
the total number of real dynamical fermion doubletsplus the hasenbusch preconditioners.

 This routine begins at "integral" time, with H and U evaluated
 at same time.
*/

#include "cl_dyn_includes.h"
#ifdef HAVE_IEEEFP_H
#include <ieeefp.h>    /* For "finite" */
#endif

void predict_next_psi(Real *oldtime,Real *newtime,Real *nexttime, int level);
int update_step(Real *oldtime,Real *newtime,Real *nexttime, double *fnorm, double *gnorm);


double update_gauge_step(float eps);
double returntrlogA;


int update()  {
  int iters=0;
  Real final_rsq;
  Real CKU0 = kappa*clov_c/(u0*u0*u0);
  /* double junktrlogA; */
  Real cg_time[2], old_cg_time[2], next_cg_time[2];
  double starttrlogA;
#ifdef HMC_ALGORITHM
  double startaction = 0, endaction, endtrlogA, change;
  Real xrandom;
#endif
  double gnorm=0.0;
  double fnorm[2];
  fnorm[0]=fnorm[1]=0.0;

#ifdef NHYP_JACOBI
  int i;
  for(i=0; i<JACOBI_HIST_MAX ; i++) jacobi_hist[i]=0 ;  
#endif  
  
#ifndef PHI_ALGORITHM
  /* quit! */
  node0_printf("only works for phi algorithm\n");
  exit(1);
#endif
  
  /* refresh the momenta */
  ranmom();
  
  
  /* higher rep code:
     DIMFxDIMF link created from NCOLxNCOL linkf after each update,
     then DIMF gauge field is switched to antiperiodic b.c.
     or to SF b.c. as needed, in fermion_rep() .
     initial DIMF links are set in update() , see also sf_make_boundary.c .
  */
  
  /* generate a pseudofermion configuration only at start*/
  make_clov(CKU0);
#ifdef LU
  starttrlogA = make_clovinv(ODD);
#else
  starttrlogA = (double)0.0;
#endif /*LU*/
  grsource_w();
  old_cg_time[0] = cg_time[0] = -1.0e6;
  old_cg_time[1] = cg_time[1] = -1.0e6;

  /* do conjugate gradient to get (Madj M)inverse * chi */
  iters += congrad_cl_m(niter,rsqmin,&final_rsq,F_OFFSET(chi[0]),F_OFFSET(psi[0]),0.0);
  if(num_masses==2)
    iters += congrad_cl_m(niter,rsqmin,&final_rsq,F_OFFSET(chi[1]),F_OFFSET(psi[1]),shift);
  cg_time[0]=cg_time[1] = 0.0;
  /**checkmul();**/
#ifdef HMC_ALGORITHM
  /* find action */
  startaction=d_action();
  startaction -= (double)2.0 * starttrlogA;
  /* printf("startaction= %g\n",startaction); */
  /* copy link field to old_link */
  gauge_field_copy_f(F_OFFSET(linkf[0]), F_OFFSET(old_linkf[0]));
#endif /*hmc*/
  
  /* do  microcanonical updating  */
  
#ifdef NHYP
  sf_coupling_flag=FORCE;
#endif
  
  iters +=  update_step(old_cg_time,cg_time,next_cg_time,
			fnorm,&gnorm);
  
#ifdef HMC_ALGORITHM
  /* find action */
  /* do conjugate gradient to get (Madj M)inverse * chi. Reuse data from update_step. */
#ifdef LU
  endtrlogA=returntrlogA;
#else
  endtrlogA = (double)0.0;
#endif 
   iters += congrad_cl_m(niter,rsqmin,&final_rsq,F_OFFSET(chi[0]),F_OFFSET(psi[0]),0.0);
  if(num_masses==2)
    iters += congrad_cl_m(niter,rsqmin,&final_rsq,F_OFFSET(chi[1]),F_OFFSET(psi[1]),shift);
#endif /* HMC */
  free_clov();  // needed for phi_alg too
#ifdef HMC_ALGORITHM
  endaction=d_action();
  endaction -= (double)2.0 * endtrlogA;
  /* printf("endaction= %g\n",endaction); */
  change = endaction-startaction;
  /* Reject configurations giving overflow */
#ifndef HAVE_IEEEFP_H
  if(fabs((double)change)>1e20){
#else
    if(!finite((double)change)){
#endif
      if(this_node==0)printf(
			     "WARNING: Correcting Apparent Overflow: Delta S = %e\n", change);
      change = 1.0e20;
    }
    
    /* decide whether to accept, if not, copy old link field back */
    /* careful - must generate only one random number for whole lattice */
    if(this_node==0)xrandom = myrand(&node_prn);
    broadcast_float(&xrandom);
    if( exp( -change ) < (double)xrandom ){
      if(traj_length > 0.0){
	gauge_field_copy_f( F_OFFSET(old_linkf[0]), F_OFFSET(linkf[0]) );
	fermion_rep();  /* BS */
      }
      if(this_node==0)printf("REJECT: delta S = %e, start S = %.12e, end S = %.12e\n",
			     change, startaction, endaction);
    }
    else {
      if(this_node==0)printf("ACCEPT: delta S = %e, start S = %.12e, end S = %.12e\n",
			     change, startaction, endaction);
    }
#endif /*HMC*/
    
    
    if(traj_length > 0){
      node0_printf("IT_PER_TRAJ %d\n", iters );
      node0_printf("MONITOR_FORCE_GAUGE %e\n",gnorm/(double)(4*nsteps[0]*nsteps[1]) );
      node0_printf("MONITOR_FORCE_FERMION0 %e\n",fnorm[0]/(double)(2*nsteps[0]) );
      node0_printf("MONITOR_FORCE_FERMION1 %e\n",fnorm[1]/(double)(4*nsteps[0]*nsteps[1]) );
#ifdef NHYP_JACOBI
      jacobi_total = 0 ;
      jacobi_avrg = 0. ;
      for(i=0; i<JACOBI_HIST_MAX ; i++){
	jacobi_total = jacobi_total + jacobi_hist[i] ;
	jacobi_avrg = jacobi_avrg + (Real)((i+1)*jacobi_hist[i]) ;
      }
      jacobi_avrg /= (Real)jacobi_total;
//    node0_printf("MONITOR_JACOBI_TOTL %d\n", jacobi_total);
      node0_printf("MONITOR_JACOBI_AVRG %.4f\nMONITOR_JACOBI_HIST   ", 
                    jacobi_avrg);
      for(i=0; i<JACOBI_HIST_MAX ; i++){
          node0_printf("%.4f ",  (Real)jacobi_hist[i]/(Real)jacobi_total );
      }
      node0_printf("\n");      
//    node0_printf("MONITOR_JACOBI_hist   ");
//    for(i=0; i<JACOBI_HIST_MAX ; i++){
//        node0_printf("%d ",  jacobi_hist[i] );
//    }
//    node0_printf("\n");
#endif
      return (iters);
    }
    
    else return(-99);
}
  
#ifdef LU
#define FORMYSITESDOMAIN FOREVENSITESDOMAIN
#else
#define FORMYSITESDOMAIN FORALLSITESDOMAIN
#endif
/* use linear extrapolation to predict next conjugate gradient solution */
/* only need even sites */
void predict_next_psi(Real *oldtime,Real *newtime,Real *nexttime, int level)
{
register int i;
register site *s;
register Real x;
wilson_vector tvec;

    if( newtime[level] != oldtime[level] ) 
      x = (nexttime[level]-newtime[level])/(newtime[level]-oldtime[level]);
    else x = 0.0;
    if( oldtime[level] < 0.0 ){
	FORMYSITESDOMAIN(i,s){
	    s->old_psi[level] = s->psi[level];
	}
    }
    else  {
	FORMYSITESDOMAIN(i,s){
            sub_wilson_vector( &(s->psi[level]), &(s->old_psi[level]), &tvec);
	    s->old_psi[level] = s->psi[level];
            scalar_mult_add_wvec( &(s->psi[level]), &tvec,x, &(s->psi[level]) );
	}
    }
    oldtime[level] = newtime[level];
    newtime[level] = nexttime[level];
}

 int update_step( Real *old_cg_time,Real *cg_time,Real *next_cg_time, double *fnorm, double *gnorm)  {
  int iters=0;
  Real CKU0 = kappa*clov_c/(u0*u0*u0);
  Real final_rsq;
  float f_eps1,f_eps0,g_eps;
  int i_multi0,i_multi1;
  Real mshift;
  int level;

  f_eps0=traj_length/(float)nsteps[0];
  f_eps1=f_eps0/2.0/(float)nsteps[1];
  g_eps=f_eps1/2.0/(float)nsteps[MAX_MASSES];
  if(num_masses==1) {level=0; mshift=0.0;}else {level=1; mshift=shift;}
  /*
  node0_printf("level %d mshift %e\n",level,mshift);
  node0_printf("f_eps0 %e f_eps1 %e g_eps %e\n",f_eps0,f_eps1,g_eps);
  */
  fnorm[0] += fermion_force(f_eps1*INT_LAMBDA,f_eps0*INT_LAMBDA);


  for(i_multi0=1;i_multi0<=nsteps[0];i_multi0++)
    {
      for(i_multi1=1;i_multi1<=nsteps[1];i_multi1++)
	{	
	  *gnorm += update_gauge_step(g_eps);
	  
	  /* save conjugate gradient solution, predict next one */
	  next_cg_time[level] = cg_time[level] + f_eps1;
	  predict_next_psi(old_cg_time,cg_time,next_cg_time,level);
	  free_clov();
	  make_clov(CKU0);
#ifdef LU
	  /* junktrlogA = */ 
	  returntrlogA=make_clovinv(ODD);
#endif /*LU*/
	  /* do conjugate gradient to get (Madj M)inverse * chi */
	  iters += congrad_cl_m(niter,rsqmin,&final_rsq,F_OFFSET(chi[level]),
				F_OFFSET(psi[level]),mshift);

	  fnorm[1] += fermion_force(f_eps1*INT_LAMBDA_MID,0.0);
	  
	  *gnorm += update_gauge_step(g_eps);

	  if(i_multi1<nsteps[1]){
	    /* save conjugate gradient solution, predict next one */
	    next_cg_time[level] = cg_time[level] + f_eps1;
	    predict_next_psi(old_cg_time,cg_time,next_cg_time,level);
	    free_clov();
	    make_clov(CKU0);
#ifdef LU
	    /* junktrlogA = */ 
	    returntrlogA=make_clovinv(ODD);
#endif /*LU*/
	    /* do conjugate gradient to get (Madj M)inverse * chi */
	    iters += congrad_cl_m(niter,rsqmin,&final_rsq,F_OFFSET(chi[level]),
				  F_OFFSET(psi[level]),mshift);
	    fnorm[1] += fermion_force(f_eps1*INT_LAMBDA_CONT,0.0);
	  }
	} /* nsteps[1] */
      
      
      /* invert both */
      
      
      next_cg_time[level] = cg_time[level] + f_eps1;
      predict_next_psi(old_cg_time,cg_time,next_cg_time,level);
      free_clov();
      make_clov(CKU0);
#ifdef LU
      /* junktrlogA = */ 
      returntrlogA=make_clovinv(ODD);
#endif /*LU*/
      /* do conjugate gradient to get (Madj M)inverse * chi */
      iters += congrad_cl_m(niter,rsqmin,&final_rsq,F_OFFSET(chi[level]),
			    F_OFFSET(psi[level]),mshift);
      if(num_masses==2){
	next_cg_time[0] = cg_time[0] + f_eps0;
	predict_next_psi(old_cg_time,cg_time,next_cg_time,0);
	iters += congrad_cl_m(niter,rsqmin,&final_rsq,F_OFFSET(chi[0]),
			      F_OFFSET(psi[0]),0.0);
      }
      
      
      fnorm[0] += fermion_force(f_eps1*INT_LAMBDA_CONT,f_eps0*INT_LAMBDA_MID);
      
      for(i_multi1=1;i_multi1<=nsteps[1];i_multi1++)
	{	
	  *gnorm += update_gauge_step(g_eps);
	  
	  /* save conjugate gradient solution, predict next one */
	  next_cg_time[level] = cg_time[level] + f_eps1;
	  predict_next_psi(old_cg_time,cg_time,next_cg_time,level);
	  /* do conjugate gradient to get (Madj M)inverse * chi */
	  free_clov();
	  make_clov(CKU0);
#ifdef LU
	  /* junktrlogA = */ 
	  returntrlogA=make_clovinv(ODD);
#endif /*LU*/
	  iters += congrad_cl_m(niter,rsqmin,&final_rsq,F_OFFSET(chi[level]),
				F_OFFSET(psi[level]),mshift);
	  fnorm[1] += fermion_force(f_eps1*INT_LAMBDA_MID,0.0);
	  
	  *gnorm += update_gauge_step(g_eps);
	  if(i_multi1<nsteps[1]){
	    /* save conjugate gradient solution, predict next one */
	    next_cg_time[level] = cg_time[level] + f_eps1;
	    predict_next_psi(old_cg_time,cg_time,next_cg_time,level);
	    /* do conjugate gradient to get (Madj M)inverse * chi */
	    free_clov();
	    make_clov(CKU0);
#ifdef LU
	    /* junktrlogA = */ 
	    returntrlogA=make_clovinv(ODD);
#endif /*LU*/
	    iters += congrad_cl_m(niter,rsqmin,&final_rsq,F_OFFSET(chi[level]),
				  F_OFFSET(psi[level]),mshift);
	    fnorm[1] += fermion_force(f_eps1*INT_LAMBDA_CONT,0.0);
	    
	  }
	}
      
            /* invert both */
      
      next_cg_time[level] = cg_time[level] + f_eps1;
      predict_next_psi(old_cg_time,cg_time,next_cg_time,level);
      /* do conjugate gradient to get (Madj M)inverse * chi */
      free_clov();
      make_clov(CKU0);
#ifdef LU
      /* junktrlogA = */ 
      returntrlogA=make_clovinv(ODD);
#endif /*LU*/
      iters += congrad_cl_m(niter,rsqmin,&final_rsq,F_OFFSET(chi[level]),
			    F_OFFSET(psi[level]),mshift);
      if(num_masses==2){
	next_cg_time[0] = cg_time[0] + f_eps0;
	predict_next_psi(old_cg_time,cg_time,next_cg_time,0);
	iters += congrad_cl_m(niter,rsqmin,&final_rsq,F_OFFSET(chi[0]),
			      F_OFFSET(psi[0]),0.0);
      }
      
      
      if(i_multi0<nsteps[0]){
	fnorm[0] += fermion_force(f_eps1*INT_LAMBDA_CONT,f_eps0*INT_LAMBDA_CONT);
      }
      else {	
	fnorm[0] += fermion_force(f_eps1*INT_LAMBDA,f_eps0*INT_LAMBDA);
      }
    } /* i_multi0 */



  return (iters);
  
} /* update_step */

/* Omelyan version;  ``dirty'' speeded-up version */
 double update_gauge_step(float eps){
  int nsw=nsteps[MAX_MASSES];
  int isw;
  double gauge_force(Real eps);
  double norm=0;
  /*
  node0_printf("gauge %d steps %e dt\n",nsw,eps);
  */
    norm+= gauge_force(eps*INT_LAMBDA);
  for(isw=1;isw<=nsw;isw++){
    update_u(0.5*eps);
    norm+= gauge_force(eps*INT_LAMBDA_MID);
    update_u(0.5*eps);
    if(isw<nsw) norm+= gauge_force(eps*INT_LAMBDA_CONT);
	else norm+= gauge_force(eps*INT_LAMBDA);
  }
  
  /* reunitarize the gauge field */
  reunitarize();
  fermion_rep();	/* BQS */
  return(norm/nsw);
}

