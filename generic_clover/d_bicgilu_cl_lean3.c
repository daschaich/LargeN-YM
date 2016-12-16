/******* d_bicgilu_cl_lean.c - BiCGstab-ILU for  clover fermions ****/
/* MIMD version 7 */

/* Modifications:
   7/18/01 calls dslash_special CD
   1/24/00 combined with Schroedinger functional version - UMH
   4/26/98 Moved parameters to structures CD
   8/29/97 ANSI prototyping and added comments C. D.
   Code created by U.M.H.
   */

/* Memory stingy version
   "r" overwrites src on even sites
   "p" overwrites src on odd sites
   3/29/00 EVENFIRST is the rule now. CD.
   */

/* Requires qic wilson vector temporaries wv1, wv2, wv3, wv4 */

/* ------------------------------------------------------------
   The matrix to be inverted is written in block even-odd form as


   M = ( R_o     -K D_oe )
       ( -K D_eo  R_e    )

   where R_o and R_e are 1 - K Clov_c/u_0^3 i sigma_mu,nu F_mu,nu
   are the site-diagonal hermitian clover matrices on odd and even
   sites.

   The dslash operators D_oe src and D_eo src are given by

   SUM_dirs ( 
      ( 1 + gamma[dir] ) * U(x,dir) * src(x+dir)
    + ( 1 - gamma[dir] ) * U_adj(x-dir,dir) * src(x-dir)
   )

   with gammas defined in libraries/mb_gamma.c.

   The ILU decomposition results in M = L A U


   M = ( 1            0 ) ( R_o  0   ) ( 1  -K/R_o D_oe )
       ( -K D_eo/R_o  1 ) ( 0    M_e ) ( 0     1        )

   where

         M_e = R_e - K^2 D_eo/R_o D_oe

   acts only on even sites.  Since M_e is not hermitian, it is
   necessary to invert M_e_dag M_e where M_e_dag = M_e^\dagger

 ------------------------------------------------------------ */

/* The source vector is in "src", and the initial guess and answer
   in "dest".  "r" is the residual vector, which is a pointer to src since
   the source  is overwritten to save space 
   and "p" and "my_mp" are
   working vectors for the conjugate gradient. 
   MinCG = minimum number of iterations.
   MaxCG = maximum number of iterations.
   size_r = desired residual, quit when we reach it.
   (Square root def for residue size_r = sqrt(r*r))
   ILU resides on parity=EVEN so do only even sites
   */

#include "generic_clover_includes.h"

/*#define CGTIME */

int bicgilu_cl(          /* Return value is number of iterations taken */
    field_offset src,    /* type wilson_vector (source vector - OVERWRITTEN!)*/
    field_offset dest,   /* type wilson_vector (answer and initial guess )*/
    quark_invert_control *qic, /* parameters controlling inversion */
    void *dmp            /* parameters defining the Dirac matrix */
    )
{
  /* Unpack required members of the structures */
  int MinCG = qic->min;      /* minimum number of iterations */
  int MaxCG = qic->max;      /* maximum number of iterations */
  Real RsdCG = qic->resid;  /* desired residual - 
				 normalized as sqrt(r*r)/sqrt(src_e*src_e */
  int flag = qic->start_flag;   /* 0: use a zero initial guess; 1: use dest */
  wilson_vector *tmp=NULL,*my_mp=NULL,*rv=NULL,*sss=NULL,*r=NULL,
    *p=NULL,*ttt=NULL,*t_dest=NULL;
  int static first_congrad = 1;

  dirac_clover_param *dcp 
    = (dirac_clover_param *)dmp; /* Cast pass-through pointer */

  Real Kappa = dcp->Kappa;     /* hopping */
  Real Clov_c = dcp->Clov_c;   /* Perturbative clover coeff */
  Real U0 = dcp->U0;           /* Tadpole correction to Clov_c */
  /* End of unpacking required members of structures */

  int N_iter;
  register int i;
  register site *s;
  Real size_src;
  double rsq, tsq;
  complex ctmp, a, b;
  complex omega, omegam;
  double_complex tdots,rvro,rvv,rvr;
  register Real MKsq = -Kappa*Kappa;
  Real CKU0 = Kappa*Clov_c/(U0*U0*U0);
#ifdef CGTIME
  double dtime;
#endif
  msg_tag *tago[8],*tage[8];
  int is_startedo, is_startede;
  
  is_startedo = is_startede = 0;

    if(even_sites_on_node!=odd_sites_on_node){
      printf("Need same number of even and odd sites on each node\n");
      terminate(1);
    }
  
  make_clov(CKU0);

  /* Take the inverse on the odd sublattice */
  make_clovinv(ODD);
  
    /* now we can allocate temporary variables and copy then */
    /* PAD may be used to avoid cache trashing */
#define PAD 0
    
    if(first_congrad) {
      tmp    = (wilson_vector *) malloc((sites_on_node+PAD)*sizeof(wilson_vector));
      my_mp  = (wilson_vector *) malloc((sites_on_node+PAD)*sizeof(wilson_vector));
      rv     = (wilson_vector *) malloc((sites_on_node+PAD)*sizeof(wilson_vector));
      sss    = (wilson_vector *) malloc((sites_on_node+PAD)*sizeof(wilson_vector));
      r      = (wilson_vector *) malloc((sites_on_node+PAD)*sizeof(wilson_vector));
      t_dest = (wilson_vector *) malloc((sites_on_node+PAD)*sizeof(wilson_vector));
      p      = r  + even_sites_on_node;
      ttt    = rv + even_sites_on_node;
      first_congrad = 0;
    }

  /* BiCGstab_ILU: */
  
  /* Start Inversion */
  
  /* src = L^(-1)*src */
#ifdef CGTIME
  dtime = -dclock();
#endif
  
  /* now we copy dest and src to temporaries */
  FORALLSITES(i,s) {
    t_dest[i] = *(wilson_vector *)F_PT(s,dest);
    r[i] = *(wilson_vector *)F_PT(s,src);
  }
  
  mult_ldu_field(r, my_mp, ODD);
  dslash_w_field_special(my_mp, my_mp, PLUS, EVEN, tage, is_startede);
  is_startede = 1;
  
  /* Normalization  */
  rsq = 0.0;
  FOREVENSITESDOMAIN(i,s) {
    scalar_mult_add_wvec( &(r[i]), &(my_mp[i]), Kappa, &(r[i]) );
    rsq += (double)magsq_wvec( &(r[i]) );
  }
  fflush(stdout);
  g_doublesum(&rsq);
  size_src = (Real)sqrt(rsq);
  
  /** if(this_node==0)printf("beginning inversion--size_src=%e\n",
      (double)size_src); **/
  
  /* r and p overwrite the source (bad for
     dynamical fermions but good for quenched calcs) */
  /* set r = src --- nothing to do */
  
  /* Initial guess */
  /* set dest = src on odd sites, whatever else you do
     (even if we restart with a nonzero solution vector, the end of the
     subroutine rebuilds the odd component from the even one. The (trivial)
     solution of the odd component of the equation is dest = src, before
     we rotate back to the basis in which  M is not checkerboard-diagonal) */
  FORODDSITESDOMAIN(i,s) {
    copy_wvec( &(r[i]), &(t_dest[i]) );
  }
  
  
  /* code if you want to start with dest=0... */
  if(flag == 0) {
    /** if(this_node==0)printf("dest_0=0\n"); **/
    FOREVENSITESDOMAIN(i,s) {
      clear_wvec( &(t_dest[i]) );
    }
  }
  /* code if you want to start dest with some particular starting value... */
  /* r=src[1]-[L^(-1)*M*U^(-1)]*dest */
  if(flag != 0) {
    /** if(this_node==0)    printf("dest_0  !=0\n"); **/
    /* we use my_mp temporarily to construct r */
    mult_ldu_field(t_dest, tmp, EVEN);
    dslash_w_field_special(t_dest, my_mp, PLUS, ODD, tago, is_startedo);
    is_startedo = 1;
    mult_ldu_field(my_mp, tmp, ODD);
    dslash_w_field_special(tmp, my_mp, PLUS, EVEN, tage, is_startede);
    is_startede = 1;
    FOREVENSITESDOMAIN(i,s) {
      scalar_mult_add_wvec( &(tmp[i]), &(my_mp[i]), MKsq, &(my_mp[i]) );
      scalar_mult_add_wvec( &(r[i]),   &(my_mp[i]), -1.0, &(r[i])     );
    }
  }
  
  rsq = 0.0;
  FOREVENSITESDOMAIN(i,s) {
    rsq += (double)magsq_wvec( &(r[i]) );
    copy_wvec( &(r[i]), &(p[i]) );
    copy_wvec( &(r[i]), &(rv[i]) );
  }
  g_doublesum(&rsq);
  qic->size_r = (Real)sqrt(rsq)/size_src;
  rvr = dcmplx(rsq,(double)0.0);
  /** if(this_node==0)    printf("beginning inversion--size_r=%e\n",
      (double)(qic->size_r)); **/
  
  for( N_iter = 0; N_iter < MinCG || (N_iter < MaxCG && RsdCG  < qic->size_r); 
      N_iter = N_iter + 1) {
    
    /*   my_mp = M(u)*p */
    mult_ldu_field(p, tmp, EVEN);
    dslash_w_field_special(p, my_mp, PLUS, ODD, tago, is_startedo);
    is_startedo = 1;
    mult_ldu_field(my_mp, tmp, ODD);
    dslash_w_field_special(tmp, my_mp, PLUS, EVEN, tage, is_startede);
    is_startede = 1;
    
    /* rvv = <rv|my_mp> */
    rvv = dcmplx((double)0.0,(double)0.0);
    FOREVENSITESDOMAIN(i,s) {
      scalar_mult_add_wvec( &(tmp[i]), &(my_mp[i]), MKsq, &(my_mp[i]) );
      ctmp = wvec_dot( &(rv[i]), &(my_mp[i]) );
      CSUM(rvv,ctmp);
    }
    g_dcomplexsum(&rvv);
    CDIV(rvr,rvv,a);
    
    /* sss = r - a*my_mp  */
    CMULREAL(a,-1.0,ctmp);
    FOREVENSITESDOMAIN(i,s) {
      c_scalar_mult_add_wvec( &(r[i]), &(my_mp[i]), &ctmp, &(sss[i]) );
    }
    
    /* ttt = M(u)*sss */
    mult_ldu_field(sss, tmp, EVEN);
    dslash_w_field_special(sss, sss, PLUS, ODD, tago, is_startedo);
    is_startedo = 1;
    mult_ldu_field(sss, tmp, ODD);
    dslash_w_field_special(tmp, ttt, PLUS, EVEN, tage, is_startede);
    is_startede = 1;
    
    /* tdots = <ttt|sss>; tsq=|ttt|^2 */
    tdots = dcmplx((double)0.0,(double)0.0);
    tsq = 0.0;
    FOREVENSITESDOMAIN(i,s) {
      scalar_mult_add_wvec( &(tmp[i]), &(ttt[i]), MKsq, &(ttt[i]) );
      ctmp = wvec_dot( &(ttt[i]), &(sss[i]) );
      CSUM(tdots, ctmp);
      tsq += (double)magsq_wvec( &(ttt[i]) );
    }
    g_dcomplexsum(&tdots);
    g_doublesum(&tsq);
    
    omega.real = tdots.real/tsq;
    omega.imag = tdots.imag/tsq;
    
    /* dest = dest + omega*sss + a*p */
    /* r = sss - omega*ttt */
    omegam.real = -omega.real;
    omegam.imag = -omega.imag;
    rsq = 0.0;
    rvro = rvr;
    rvr = dcmplx((double)0.0,(double)0.0);
    FOREVENSITESDOMAIN(i,s) {
      c_scalar_mult_add_wvec( &(t_dest[i]), &(sss[i]), &omega,  &(t_dest[i]) );
      c_scalar_mult_add_wvec( &(t_dest[i]), &(p[i]),   &a,      &(t_dest[i]) );
      c_scalar_mult_add_wvec( &(sss[i]),    &(ttt[i]), &omegam, &(r[i])      );
      ctmp = wvec_dot( &(rv[i]), &(r[i]) );
      CSUM(rvr, ctmp);
      rsq += (double)magsq_wvec( &(r[i]) );
    }
    g_dcomplexsum(&rvr);
    g_doublesum(&rsq);
    
    CDIV(rvr,rvro,b);
    CMUL(a,b,ctmp);
    CDIV(ctmp,omega,b);
    
    /*   p = r + b*(p - omega*my_mp)  */
    FOREVENSITESDOMAIN(i,s) {
      c_scalar_mult_add_wvec( &(p[i]),  &(my_mp[i]), &omegam, &(p[i]) );
      c_scalar_mult_add_wvec( &(r[i]),  &(p[i]),     &b,      &(p[i]) );
    }
    
    qic->size_r = (Real)sqrt(rsq)/size_src;
    /** if(this_node==0){printf("iteration= %d, residue= %e\n",N_iter,
	(double)(qic->size_r));fflush(stdout); **/
  }
#ifdef CGTIME
  dtime += dclock();
#endif
  if(this_node==0){
    if(N_iter==0)
      printf("BiCGILU: NO iterations taken size_r= %.2e\n",qic->size_r);
#ifdef CGTIME
    else
      printf("BiCGILU: time = %.2e size_r= %.2e iters= %d MF = %.1f\n",
	     dtime,qic->size_r,N_iter,
	     (double)8742*N_iter*even_sites_on_node/(dtime*1e6));
#endif
    fflush(stdout);
  }
  
  /** if( (qic->size_r) > RsdCG ) {
    if(this_node==0){printf(" BiCG_ILU: Not Converged\n");fflush(stdout);}
    } **/
  
  /* dest = R^(-1)*dest  */
  dslash_w_field_special(t_dest, my_mp, PLUS, ODD, tago, is_startedo);
  is_startedo = 1;
  FORODDSITESDOMAIN(i,s) {
    scalar_mult_add_wvec( &(t_dest[i]), &(my_mp[i]), Kappa, &(my_mp[i]) );
  }
  mult_ldu_field(my_mp, t_dest, ODD);
  
  /* Copy back for all sites and free memory */
  FORALLSITES(i,s) {
    *(wilson_vector *)F_PT(s,dest) = t_dest[i];
  }
  free(tmp); free(my_mp); free(rv); free(sss); free(r); free(t_dest);
  first_congrad = 1;

  for( i=XUP; i <= TUP; i++) {
    if(is_startedo)cleanup_gather(tago[i]);
    if(is_startedo)cleanup_gather(tago[OPP_DIR(i)]);
    if(is_startede)cleanup_gather(tage[i]);
    if(is_startede)cleanup_gather(tage[OPP_DIR(i)]);
  }
  is_startede = is_startedo = 0;
    
  free_clov();
  cleanup_tmp_links();
  return(N_iter);
}

