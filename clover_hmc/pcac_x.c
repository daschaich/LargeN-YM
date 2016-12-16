/*
* wilson/clover fermions

* Fermionic observables intended for screening masses under SF.
  Can be used *both* with or without SF.
  To be called from control.c provided nx>nt.

* clover operator:
* M = A - kappa*( Dslash_eo + DSLASH_oe )
* 2*kappa=1/(4+m0)

* OBSERVABLE                                   FLAGS (if any)

XXPPN = -<pi(x) pi(0)>   "normal" b.c.
XXAPN = <A_1(x) pi(0)>   "normal" b.c.
XXVVN = <V_1(x) V_1(0)>  "normal" b.c.         V_TOO
XXPPD = -<pi(x) pi(0)>   "compouand" b.c.      D_TOO
XXAPD = <A_1(x) pi(0)>   "compouand" b.c.      D_TOO
XXVVD = <V_1(x) V_1(0)>  "compouand" b.c.      D_TOO, V_TOO

COMPILATION FLAGS controlled by defines.h
Measurements at py=pz=0.
Nominally also pt=0, but (possibly) with SF boundary conditions.

"normal" b.c. =
CL: periodic
SF: ferm_twist_phase

"compouand" b.c. =
CL: linear combination of p.b.c. and a.p.b.c.
SF: linear combination of ferm_twist_phase and ferm_twist_phase+Pi

See long comment below for detailed form of observables.

***********************************************************************/

#include "cl_dyn_includes.h"

/* Assumes LU preconditioned in congrad */
#ifndef LU
BOMB THE COMPILE
#endif

int pcac_x() {

   register int i,j,k,xx;
   register site *s;
   register complex cc;
   int iters=0;
#ifdef SF
   int yztvolume=(nt-1)*ny*nz;
#else
   int yztvolume=nt*ny*nz;
#endif
   complex *sinkpn,*sinkan;   /* nx dot products         */
   Real *corrppn,*corrapn;  /* nx results                */
#ifdef D_TOO
   complex *sinkpd,*sinkad;   /* nx dot products         */
   Real *corrppd,*corrapd;  /* nx results                */
#endif
#ifdef V_TOO
   complex *sinkvn,*srcvn;   /* nx dot products          */
   Real *corrvvn;  /* nx results                         */
#ifdef D_TOO
   complex *sinkvd,*srcvd;   /* nx dot products          */
   Real *corrvvd;  /* nx results                         */
#endif
#endif

   sinkpn = (complex *)malloc(nx*sizeof(complex));
   sinkan = (complex *)malloc(nx*sizeof(complex));
   corrppn = (Real *)malloc(nx*sizeof(Real));
   corrapn = (Real *)malloc(nx*sizeof(Real));
#ifdef D_TOO
   sinkpd = (complex *)malloc(nx*sizeof(complex));
   sinkad = (complex *)malloc(nx*sizeof(complex));
   corrppd = (Real *)malloc(nx*sizeof(Real));
   corrapd = (Real *)malloc(nx*sizeof(Real));
#endif
#ifdef V_TOO
   sinkvn = (complex *)malloc(nx*sizeof(complex));
   srcvn  = (complex *)malloc(nx*sizeof(complex));
   corrvvn = (Real *)malloc(nx*sizeof(Real));
#ifdef D_TOO
   sinkvd = (complex *)malloc(nx*sizeof(complex));
   srcvd  = (complex *)malloc(nx*sizeof(complex));
   corrvvd = (Real *)malloc(nx*sizeof(Real));
#endif
#endif

/***  generate gaussian random source on x-slice x=0 only  *************
*  normalization of random sources:        SEE libraries/gaussrand.c   *
*  <eta|eta> = <Re eta|Re eta> + <Im eta|Im eta> = 1/2 + 1/2 = 1      */

   FORALLSITES(i,s){
      for(k=0;k<4;k++)for(j=0;j<DIMF;j++){
#ifdef SF
         if(s->x == 0 && s->t > 0){
#else
         if(s->x == 0){
#endif
#ifdef SITERAND
            s->g_rand.d[k].c[j].real = gaussian_rand_no(&(s->site_prn));
            s->g_rand.d[k].c[j].imag = gaussian_rand_no(&(s->site_prn));
#else
            s->g_rand.d[k].c[j].real = gaussian_rand_no(&node_prn);
            s->g_rand.d[k].c[j].imag = gaussian_rand_no(&node_prn);
#endif
         }
         else{
            s->g_rand.d[k].c[j].real = 0;
            s->g_rand.d[k].c[j].imag = 0;
         }
      }
   }

/***  inversions  ******************************************************
*  "Compound b.c.: sum propagators with a.p.b.c. and with p.b.c.      */

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

#ifdef D_TOO
/* APBC inversions. Copy g_rand to avoid overwriting it       */

   boundary_flip_x(MINUS);

/* compute 1/M_apbc |src>                                   */

   FORALLSITESDOMAIN(i,s) {
      copy_wvec( &(s->g_rand), &(s->chi[0]) );
   }
   qic.start_flag = 0;   /* Use zero initial guess for dest */
   iters += wilson_invert(F_OFFSET(chi[0]),F_OFFSET(psi[0]),F_OFFSET(r),
                         CURINVERT,&qic,(void *)&dcp);

#ifdef V_TOO
/* compute 1/M_apbc |gamma_5*gamma_1*src>                   */

   FORALLSITES(i,s) {
      mult_by_gamma( &(s->g_rand), &(s->chi[0]), XUP );
      mult_by_gamma( &(s->chi[0]), &(s->chi[1]), GAMMAFIVE );
   }
   qic.start_flag = 0;   /* Use zero initial guess for dest */
   iters += wilson_invert(F_OFFSET(chi[1]),F_OFFSET(psi[1]),F_OFFSET(r),
                         CURINVERT,&qic,(void *)&dcp);
#endif /* V_TOO */

   boundary_flip_x(PLUS);
#endif /* D_TOO */

/* PBC inversions, can overwrite g_rand on last inversion   */

/* compute |invp> = 1/M_pbc |src>                           */

   FORALLSITESDOMAIN(i,s) {
      copy_wvec( &(s->g_rand), &(s->chi[0]) );
   }
   qic.start_flag = 0;   /* Use zero initial guess for dest */
   iters += wilson_invert(F_OFFSET(chi[0]),F_OFFSET(invp),F_OFFSET(r),
                         CURINVERT,&qic,(void *)&dcp);

#ifdef V_TOO
/* compute |invvp> = 1/M_pbc |gamma_5*gamma_1*src>          */

   FORALLSITES(i,s) {
      mult_by_gamma( &(s->g_rand), &(s->chi[0]), XUP );
      mult_by_gamma( &(s->chi[0]), &(s->g_rand), GAMMAFIVE );
   }
   qic.start_flag = 0;   /* Use zero initial guess for dest */
   iters += wilson_invert(F_OFFSET(g_rand),F_OFFSET(chi[1]),F_OFFSET(r),
                         CURINVERT,&qic,(void *)&dcp);
#endif /* V_TOO */

#ifdef D_TOO
/* compute 1/M_cbc |src>                       */

   FORALLSITESDOMAIN(i,s) {
      add_wilson_vector( &(s->psi[0]), &(s->invp), &(s->inva) );
      add_wilson_vector( &(s->psi[1]), &(s->chi[1]), &(s->psi[1]) );
   }
#endif /* D_TOO */

/***  observables  *****************************************************

all observables at py=pz=pt=0.
yztvolume = 3-volume
suitably modified for SF.

|src>  = random gaussian vector on x-slice = 0
|sink> = random gaussian vector on x-slice = x
|invp> = 1/M_pbc |src>
|inva> = 1/M_cbc |src>
|invvp> = 1/M_pbc |gamma_5*gamma_1*src>
|invva> = 1/M_cbc |gamma_5*gamma_1*src>

to save space:
        chi[1] functions as invvp
        psi[1] functions as invva

P ~ psibar gamma_5 psi
A_1 ~ psibar gamma_5 gamma_1 psi
V_1 ~ psibar gamma_1 psi

XXPPN = -<pi(x) pi(0)>  = <sink|invp> <invp|sink> / yztvolume
XXAPN = <A_1(x) pi(0)>  = -Re <sink*gamma_1|invp> <invp|sink> / yztvolume
XXVVN = <V_1(t) V_1(0)> = Re <sink|invvp> <invp|gamma_5*gamma_1*sink>/yztvolume
XX**D = like XX**N, except with "compound" b.c.

Output format:

PROP TOTAL_ITER iters
PROP XXPPN result(x=0) ... result(x=nx-1)
PROP XXAPN result(x=0) ... result(x=nx-1)
PROP XXVVN result(x=0) ... result(x=nx-1)
PROP XXPPD result(x=0) ... result(x=nx-1)
PROP XXAPD result(x=0) ... result(x=nx-1)
PROP XXVVD result(x=0) ... result(x=nx-1)

*/

/***  sink gaussian random vector  ************************************/

   FORALLSITESDOMAIN(i,s){
      for(k=0;k<4;k++)for(j=0;j<DIMF;j++){
#ifdef SITERAND
            s->g_rand.d[k].c[j].real = gaussian_rand_no(&(s->site_prn));
            s->g_rand.d[k].c[j].imag = gaussian_rand_no(&(s->site_prn));
#else
            s->g_rand.d[k].c[j].real = gaussian_rand_no(&node_prn);
            s->g_rand.d[k].c[j].imag = gaussian_rand_no(&node_prn);
#endif
      }
   }

/***  accummulate dot products  ***************************************/

   for(xx=0;xx<nx;xx++){
      sinkpn[xx] = cmplx(0.0,0.0);
      sinkan[xx] = cmplx(0.0,0.0);
#ifdef D_TOO
      sinkpd[xx] = cmplx(0.0,0.0);
      sinkad[xx] = cmplx(0.0,0.0);
#endif
#ifdef V_TOO
      sinkvn[xx] = cmplx(0.0,0.0);
      srcvn[xx] = cmplx(0.0,0.0);
#ifdef D_TOO
      sinkvd[xx] = cmplx(0.0,0.0);
      srcvd[xx] = cmplx(0.0,0.0);
#endif
#endif
   }

   FORALLSITESDOMAIN(i,s){
      cc = wvec_dot( &(s->g_rand), &(s->invp) );
      CSUM(sinkpn[s->x],cc);
      mult_by_gamma( &(s->g_rand), &(s->chi[0]), XUP );
      cc = wvec_dot( &(s->chi[0]), &(s->invp) );
      CSUM(sinkan[s->x],cc);
#ifdef D_TOO
      cc = wvec_dot( &(s->g_rand), &(s->inva) );
      CSUM(sinkpd[s->x],cc);
      cc = wvec_dot( &(s->chi[0]), &(s->inva) );
      CSUM(sinkad[s->x],cc);
#endif /* D_TOO */
#ifdef V_TOO
      cc = wvec_dot( &(s->g_rand), &(s->chi[1]) );
      CSUM(sinkvn[s->x],cc);
      mult_by_gamma( &(s->g_rand), &(s->psi[0]), XUP );
      mult_by_gamma( &(s->psi[0]), &(s->chi[0]), GAMMAFIVE );
      cc = wvec_dot( &(s->chi[0]), &(s->invp) );
      CSUM(srcvn[s->x],cc);
#ifdef D_TOO
      cc = wvec_dot( &(s->g_rand), &(s->psi[1]) );
      CSUM(sinkvd[s->x],cc);
      mult_by_gamma( &(s->g_rand), &(s->psi[0]), XUP );
      mult_by_gamma( &(s->psi[0]), &(s->chi[0]), GAMMAFIVE );
      cc = wvec_dot( &(s->chi[0]), &(s->inva) );
      CSUM(srcvd[s->x],cc);
#endif /* D_TOO */
#endif /* V_TOO */
   }

   for(xx=0;xx<nx;xx++){
      g_complexsum( &(sinkpn[xx]) );
      g_complexsum( &(sinkan[xx]) );
#ifdef D_TOO
      g_complexsum( &(sinkpd[xx]) );
      g_complexsum( &(sinkad[xx]) );
#endif /* D_TOO */
#ifdef V_TOO
      g_complexsum( &(sinkvn[xx]) );
      g_complexsum( &(srcvn[xx]) );
#ifdef D_TOO
      g_complexsum( &(sinkvd[xx]) );
      g_complexsum( &(srcvd[xx]) );
#endif /* D_TOO */
#endif /* V_TOO */
   }

/***  compute correlators  ********************************************/

   for(xx=0;xx<nx;xx++){
      corrppn[xx] = sinkpn[xx].real*sinkpn[xx].real
                   +sinkpn[xx].imag*sinkpn[xx].imag;
      corrppn[xx] *= (4*kappa*kappa)/yztvolume;

      corrapn[xx] = sinkan[xx].real*sinkpn[xx].real
                   +sinkan[xx].imag*sinkpn[xx].imag;
      corrapn[xx] *= (-4*kappa*kappa)/yztvolume;
#ifdef D_TOO
      corrppd[xx] = sinkpd[xx].real*sinkpd[xx].real
                   +sinkpd[xx].imag*sinkpd[xx].imag;
      corrppd[xx] *= (kappa*kappa)/yztvolume;

      corrapd[xx] = sinkad[xx].real*sinkpd[xx].real
                   +sinkad[xx].imag*sinkpd[xx].imag;
      corrapd[xx] *= (-kappa*kappa)/yztvolume;
#endif /* D_TOO */
#ifdef V_TOO
      corrvvn[xx] = sinkvn[xx].real*srcvn[xx].real
                   +sinkvn[xx].imag*srcvn[xx].imag;
      corrvvn[xx] *= (-4*kappa*kappa)/yztvolume;
#ifdef D_TOO
      corrvvd[xx] = sinkvd[xx].real*srcvd[xx].real
                   +sinkvd[xx].imag*srcvd[xx].imag;
      corrvvd[xx] *= (-4*kappa*kappa)/yztvolume;
#endif /* D_TOO */
#endif /* V_TOO */
   }

/***  printouts  ******************************************************/

   if(this_node==0) {
      printf("PROP TOTAL_ITERS %d\n", iters);

/* print x-slices correlations                */

      printf("PROP XXPPN");
      for(xx=0;xx<nx;xx++){
         printf(" %e", (double)corrppn[xx]);
      }
      printf("\n");

      printf("PROP XXAPN");
      for(xx=0;xx<nx;xx++){
         printf(" %e", (double)corrapn[xx]);
      }
      printf("\n");

#ifdef V_TOO
      printf("PROP XXVVN");
      for(xx=0;xx<nx;xx++){
         printf(" %e", (double)corrvvn[xx]);
      }
      printf("\n");
#endif /* V_TOO */

#ifdef D_TOO
      printf("PROP XXPPD");
      for(xx=0;xx<nx;xx++){
         printf(" %e", (double)corrppd[xx]);
      }
      printf("\n");

      printf("PROP XXAPD");
      for(xx=0;xx<nx;xx++){
         printf(" %e", (double)corrapd[xx]);
      }
      printf("\n");

#ifdef V_TOO
      printf("PROP XXVVD");
      for(xx=0;xx<nx;xx++){
         printf(" %e", (double)corrvvd[xx]);
      }
      printf("\n");
#endif /* V_TOO */
#endif /* D_TOO */
   }

   free(sinkpn); free(sinkan); free(corrppn); free(corrapn);
#ifdef D_TOO
   free(sinkpd); free(sinkad); free(corrppd); free(corrapd);
#endif /* D_TOO */
#ifdef V_TOO
   free(sinkvn); free(srcvn); free(corrvvn);
#ifdef D_TOO
   free(sinkvd); free(srcvd); free(corrvvd);
#endif /* D_TOO */
#endif /* V_TOO */

   fflush(stdout);
   return(iters);
}
/* end pcac_x */

