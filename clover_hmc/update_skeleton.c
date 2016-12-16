
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

  fnorm[0] += fermion_force(f_eps1*INT_LAMBDA,f_eps0*INT_LAMBDA);

  for(i_multi0=1;i_multi0<=nsteps[0];i_multi0++)
    {
      for(i_multi1=1;i_multi1<=nsteps[1];i_multi1++)
	{
	  *gnorm += update_gauge_step(g_eps);
	  iters += congrad_cl_m(niter,rsqmin,&final_rsq,F_OFFSET(chi[level]),
				F_OFFSET(psi[level]),mshift);
	  fnorm[1] += fermion_force(f_eps1*INT_LAMBDA_MID,0.0);
	  *gnorm += update_gauge_step(g_eps);

	  if(i_multi1<nsteps[1]){
	    iters += congrad_cl_m(niter,rsqmin,&final_rsq,F_OFFSET(chi[level]),
				  F_OFFSET(psi[level]),mshift);
	    fnorm[1] += fermion_force(f_eps1*INT_LAMBDA_CONT,0.0);
	  }
	} /* nsteps[1] */

      /* invert both */
      iters += congrad_cl_m(niter,rsqmin,&final_rsq,F_OFFSET(chi[level]),
			    F_OFFSET(psi[level]),mshift);
      if(num_masses==2){
	iters += congrad_cl_m(niter,rsqmin,&final_rsq,F_OFFSET(chi[0]),
			      F_OFFSET(psi[0]),0.0);
      }
      fnorm[0] += fermion_force(f_eps1*INT_LAMBDA_CONT,f_eps0*INT_LAMBDA_MID);

      for(i_multi1=1;i_multi1<=nsteps[1];i_multi1++)
	{
	  *gnorm += update_gauge_step(g_eps);
	  iters += congrad_cl_m(niter,rsqmin,&final_rsq,F_OFFSET(chi[level]),
				F_OFFSET(psi[level]),mshift);
	  fnorm[1] += fermion_force(f_eps1*INT_LAMBDA_MID,0.0);
	  *gnorm += update_gauge_step(g_eps);
	  if(i_multi1<nsteps[1]){
	    iters += congrad_cl_m(niter,rsqmin,&final_rsq,F_OFFSET(chi[level]),
				  F_OFFSET(psi[level]),mshift);
	    fnorm[1] += fermion_force(f_eps1*INT_LAMBDA_CONT,0.0);
	  }
	}

            /* invert both */
      iters += congrad_cl_m(niter,rsqmin,&final_rsq,F_OFFSET(chi[level]),
			    F_OFFSET(psi[level]),mshift);
      if(num_masses==2){
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
}

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

