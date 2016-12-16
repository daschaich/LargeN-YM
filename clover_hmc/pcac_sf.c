/************************** pcac_sf.c **********************************
SF fermion observables relevant for AWI and PCAC .

NOTE! hopping term in MILC conventions is opposite from alpha collaboration:
  MILC:  psibar(t) P_+ psi(t+1)
  alpha: psibar(t) P_- psi(t+1)

where P_+- = (1/2)(1 +- gamma_0).  Thus we must replace P_+ <--> P_-
in alpha's formulae before implementing in code.
See libraries/wp_shrink.c for gamma matrix conventions.

Measurements (sign conventions as in pcac_t.c):

  TTPDN = -<pi(t)  pi_sf(1)>
  TTADN =  <A_4(t) pi_sf(1)>

if UP_TOO is set in defines.h we also do:

  TTPUP = -<pi(t)  pi_sf(nt-1)>
  TTAUP = -<A_4(t) pi_sf(nt-1)>

*************************************************************************/

#include "cl_dyn_includes.h"

/* Assumes LU preconditioned in congrad */
#ifndef LU
BOMB THE COMPILE
#endif

void do_measure( double *corrp, double *corra );
void wall_to_wall();
double w2w[8];

int pcac_sf() {

  register int i,j,k,tt;
  register site *s;
  register complex cc;
  int iters=0;
/* int tmin=1; adjust as wanted */
  int svolume=nx*ny*nz;
  msg_tag *mtag;
  double *corrpdn,*corradn;   /* nt results              */
#ifdef UP_TOO
  Real *corrpup,*corraup;
#endif

  for(i=0;i<8;i++) w2w[i]=(double)0.0;

/*
For measurements only at t=tmin+1, tmin+2, ..., nt-tmin-1.
Require nt >= 2*(tmin+1).

  if( 2*(tmin+1)>nt ) {
    if(this_node==0) {
      printf("sorry, tmin=%d is larger than nt/2-1\n", tmin);
    }
    terminate(1);
  }
*/
  corrpdn = (double *)malloc((nt-1)*sizeof(double));
  corradn = (double *)malloc((nt-1)*sizeof(double));
#ifdef UP_TOO
  corrpup = (double *)malloc((nt-1)*sizeof(double));
  corraup = (double *)malloc((nt-1)*sizeof(double));
#endif

  for(tt=0;tt<nt-1;tt++){
    corrpdn[tt] = (double)0.0;
    corradn[tt] = (double)0.0;
#ifdef UP_TOO
    corrpup[tt] = (double)0.0;
    corraup[tt] = (double)0.0;
#endif
  }

#ifdef BI
#define CURINVERT bicgilu_cl
  /* Load temporaries specific to inverter */
  qic.wv1 = F_OFFSET(tmp);
  qic.wv2 = F_OFFSET(mp);
  qic.wv3 = F_OFFSET(p);
  qic.wv4 = F_OFFSET(sss);
#else
#define CURINVERT cgilu_cl
  /* Load temporaries specific to inverter */
  qic.wv1 = F_OFFSET(tmp);
  qic.wv2 = F_OFFSET(mp);
#endif

  qic.start_flag = 0;   /* Use zero initial guess for dest */

/* observables require one inversion per color per eigenvector of P_+- */

  for(j=0;j<DIMF;j++) {

/*****  NOTE!!!: gather must be redone each time,                   ***/
/*****  because each inversion "corrupts" the contents of gen_pt[0] ***/

    /*  t=1 sources; h^-_i, i=1,2, are ev's of P_-  */
    /*  |src(t=1)> = U^dag_4(t=0)|h^-_i,color=j>    */

    /*  prepare source for h^-_1 = (1, 0, -1, 0)    */
    /*  first get link[TUP] from direction TDOWN    */
    mtag = start_gather_site( F_OFFSET(link[TUP]), sizeof(su3_matrix),
                              TDOWN, EVENANDODD, gen_pt[0] );
    wait_gather(mtag);

    FORALLSITES(i,s) {
      if(s->t==1) {
        for(k=0;k<DIMF;k++) {
          CONJG( ((su3_matrix *)(gen_pt[0][i]))->e[j][k] , cc );
          s->psi[0].d[0].c[k] = cc;
          CNEGATE( cc, s->psi[0].d[2].c[k]);
          s->psi[0].d[1].c[k] = cmplx(0.0,0.0);
          s->psi[0].d[3].c[k] = cmplx(0.0,0.0);
        }
      }
      else {
        clear_wvec( &(s->psi[0]) );
      }
    }

    /* we cleanup before calling the inversion routine */
    cleanup_gather(mtag);

    iters += wilson_invert(F_OFFSET(psi[0]),F_OFFSET(invp),F_OFFSET(r),
                           CURINVERT,&qic,(void *)&dcp);
    do_measure( corrpdn, corradn );
    wall_to_wall();

    /*  prepare source for h^-_2 = (0, 1, 0, -1)    */
    /*  first get link[TUP] from direction TDOWN    */
    mtag = start_gather_site( F_OFFSET(link[TUP]), sizeof(su3_matrix),
                              TDOWN, EVENANDODD, gen_pt[0] );
    wait_gather(mtag);

    FORALLSITES(i,s) {
      if(s->t==1) {
        for(k=0;k<DIMF;k++) {
          CONJG( ((su3_matrix *)(gen_pt[0][i]))->e[j][k] , cc );
          s->psi[0].d[1].c[k] = cc;
          CNEGATE( cc, s->psi[0].d[3].c[k]);
          s->psi[0].d[0].c[k] = cmplx(0.0,0.0);
          s->psi[0].d[2].c[k] = cmplx(0.0,0.0);
        }
      }
      else {
        clear_wvec( &(s->psi[0]) );
      }
    }

    /* we cleanup before calling the inversion routine */
    cleanup_gather(mtag);

    iters += wilson_invert(F_OFFSET(psi[0]),F_OFFSET(invp),F_OFFSET(r),
                           CURINVERT,&qic,(void *)&dcp);
    do_measure( corrpdn, corradn );
    wall_to_wall();

#ifdef UP_TOO
    /*  t=nt-1 sources; h^+_i, i=1,2, are ev's of P_+  */
    /*  |src(t=nt-1)> = U_4(t=nt-1)|h^+_i,color=j>     */

    /*  prepare source for h^+_1 = (1, 0, 1, 0)        */
    FORALLSITES(i,s) {
      if(s->t==nt-1) {
        for(k=0;k<DIMF;k++) {
          cc = s->link[TUP].e[k][j];
          s->psi[0].d[0].c[k] = cc;
          s->psi[0].d[2].c[k] = cc;
          s->psi[0].d[1].c[k] = cmplx(0.0,0.0);
          s->psi[0].d[3].c[k] = cmplx(0.0,0.0);
        }
      }
      else {
        clear_wvec( &(s->psi[0]) );
      }
    }

    iters += wilson_invert(F_OFFSET(psi[0]),F_OFFSET(invp),F_OFFSET(r),
                           CURINVERT,&qic,(void *)&dcp);
    do_measure( corrpup, corraup );

    /*  prepare source for h^+_2 = (0, 1, 0, 1)        */
    FORALLSITES(i,s) {
      if(s->t==nt-1) {
        for(k=0;k<DIMF;k++) {
          cc = s->link[TUP].e[k][j];
          s->psi[0].d[1].c[k] = cc;
          s->psi[0].d[3].c[k] = cc;
          s->psi[0].d[0].c[k] = cmplx(0.0,0.0);
          s->psi[0].d[2].c[k] = cmplx(0.0,0.0);
        }
      }
      else {
        clear_wvec( &(s->psi[0]) );
      }
    }

    iters += wilson_invert(F_OFFSET(psi[0]),F_OFFSET(invp),F_OFFSET(r),
                           CURINVERT,&qic,(void *)&dcp);
    do_measure( corrpup, corraup );
#endif /* UP_TOO */

  } /* end loop over source's color */

  g_vecdoublesum( (double *)corrpdn, nt-1 );
  g_vecdoublesum( (double *)corradn, nt-1 );
#ifdef UP_TOO
  g_vecdoublesum( (double *)corrpup, nt-1 );
  g_vecdoublesum( (double *)corraup, nt-1 );
#endif

    /* normalization and sign */
  for(tt=0;tt<nt-1;tt++){
    corrpdn[tt] *= (4*kappa*kappa)/svolume;
    corradn[tt] *= (-4*kappa*kappa)/svolume;
#ifdef UP_TOO
    corrpup[tt] *= (4*kappa*kappa)/svolume;
    corraup[tt] *= (4*kappa*kappa)/svolume;
#endif
  }

  for(i=0;i<8;i++){
     w2w[i] *= 4*kappa*kappa/(svolume*svolume);
  }

  /* print out results */

   if(this_node==0) {
      printf("PROP TOTAL_ITERS %d\n", iters);

  /* printout momenta in this order:
     (000), (001), (010), (100), (011), (101), (110), (111)        */
      printf("PROP W2W %e %e %e %e %e %e %e %e\n",
             w2w[0],w2w[1],w2w[2],w2w[4],w2w[3],w2w[5],w2w[6],w2w[7]);

  /* print time-slices correlations             */

      printf("PROP PDN");
      for(tt=0;tt<nt-1;tt++){
/*       printf(" %e", (double)corrpdn[tt]); */
         printf(" %e", corrpdn[tt]);
      }
      printf("\n");

      printf("PROP ADN");
      for(tt=0;tt<nt-1;tt++){
         printf(" %e", corradn[tt]);
      }
      printf("\n");

#ifdef UP_TOO
      printf("PROP PUP");
      for(tt=0;tt<nt-1;tt++){
         printf(" %e", corrpup[tt]);
      }
      printf("\n");

      printf("PROP AUP");
      for(tt=0;tt<nt-1;tt++){
         printf(" %e", corraup[tt]);
      }
      printf("\n");
#endif

  }

  free(corrpdn);  free(corradn);
#ifdef UP_TOO
  free(corrpup);  free(corraup);
#endif

  fflush(stdout);
  return(iters);

} /* pcac_sf */


void do_measure( double *corrp, double *corra ) {
  register int i;
  register site *s;
  wilson_vector wtmp;

  FORALLSITESDOMAIN(i,s){
    corrp[s->t - 1] += (double)magsq_wvec( &(s->invp) );
    mult_by_gamma( &(s->invp), &wtmp, TUP );
    corra[s->t - 1] += (double)wvec_rdot( &(s->invp), &wtmp );
  }
} /* do_measure */

/*
Wall-to-wall correlator.
For Z_p, the following is summed over the t=1 sources:

w = sum_vec{x} U^dag_4(vec{x},t=nt-1) invp(vec{x},t=nt-1)
w2w = w^dag (1-gamma_0) w

Normalization of w2w: we use (1-gamma_0), and do not divide by 2,
to keep the same normalization as for the eigenfunctions we use
in the source.

w2w[mm]: we project at the sink on 8 different momenta,
where each momentum component px,py,pz, can be either 0 or 2*PI/L.
*/

void wall_to_wall() {
  register int i,j,mm,mmx,mmy,mmz;
  register site *s;
  Real theta;
  complex phase[8];
  wilson_vector wtmp;
  dhalf_wilson_vector hw[8];

  for(mm=0;mm<8;mm++)for(j=0;j<DIMF;j++){
    hw[mm].h[0].c[j].real = hw[mm].h[0].c[j].imag = (double)0.0;
    hw[mm].h[1].c[j].real = hw[mm].h[1].c[j].imag = (double)0.0;
  }

  FORALLSITES(i,s) {
    if(s->t==nt-1) {
/* parallel transport to t=nt boundary */
      mult_adj_mat_wilson_vec( &(s->link[TUP]), &(s->invp), &wtmp);
/* multiply by plane wave and add up,
   picking only eigenvectors of P_-    */
      for(mmx=0;mmx<2;mmx++)for(mmy=0;mmy<2;mmy++)for(mmz=0;mmz<2;mmz++){
         theta = 2*PI*( (s->x)*mmx/(Real)nx
                       +(s->y)*mmy/(Real)ny
                       +(s->z)*mmz/(Real)nz );
         phase[4*mmx+2*mmy+mmz].real = cos(theta);
         phase[4*mmx+2*mmy+mmz].imag = sin(theta);
      }
      for(mm=0;mm<8;mm++)for(j=0;j<DIMF;j++){
         hw[mm].h[0].c[j].real += (double)(phase[mm].real*wtmp.d[0].c[j].real
                                          -phase[mm].imag*wtmp.d[0].c[j].imag
                                          -phase[mm].real*wtmp.d[2].c[j].real
                                          +phase[mm].imag*wtmp.d[2].c[j].imag);
         hw[mm].h[0].c[j].imag += (double)(phase[mm].real*wtmp.d[0].c[j].imag
                                          +phase[mm].imag*wtmp.d[0].c[j].real
                                          -phase[mm].real*wtmp.d[2].c[j].imag
                                          -phase[mm].imag*wtmp.d[2].c[j].real);
         hw[mm].h[1].c[j].real += (double)(phase[mm].real*wtmp.d[1].c[j].real
                                          -phase[mm].imag*wtmp.d[1].c[j].imag
                                          -phase[mm].real*wtmp.d[3].c[j].real
                                          +phase[mm].imag*wtmp.d[3].c[j].imag);
         hw[mm].h[1].c[j].imag += (double)(phase[mm].real*wtmp.d[1].c[j].imag
                                          +phase[mm].imag*wtmp.d[1].c[j].real
                                          -phase[mm].real*wtmp.d[3].c[j].imag
                                          -phase[mm].imag*wtmp.d[3].c[j].real);
      }
    }
  }

/* also ok:  g_vecdcomplexsum( (double_complex *)hw, 16*DIMF ); */
  g_vecdoublesum( (double *)hw, 32*DIMF );

  for(mm=0;mm<8;mm++)for(j=0;j<DIMF;j++){
    w2w[mm] += hw[mm].h[0].c[j].real*hw[mm].h[0].c[j].real
              +hw[mm].h[0].c[j].imag*hw[mm].h[0].c[j].imag
              +hw[mm].h[1].c[j].real*hw[mm].h[1].c[j].real
              +hw[mm].h[1].c[j].imag*hw[mm].h[1].c[j].imag ;
  }

} /* wall_to_wall */
