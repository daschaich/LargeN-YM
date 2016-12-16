/*
* wilson/clover fermions

* Fermionic observables to find kappa_c by imposing PCAC relation:

* 4 = number of Dirac components
* NCOL = number of colors
* DIMF = size of fermion rep

* clover operator:
* M = A - kappa*( Dslash_eo + DSLASH_oe )
* 2*kappa=1/(4+m0)

* measurements (py=pz=pt=0):
* XXPP = -<pi(x) pi(0)>
* XXAP = <A_1(x) pi(0)>
* XXPA = -<pi(x) A_1(0)>
* XXAA = <A_1(x) A_1(0)>

if V_TOO is set in defines.h we also do:
* XXVV = <V_3(x) V_3(0)>

* specific members of site structure:
        wilson_vector invp;
        wilson_vector inva;

* "compouand b.c." = linear combination of p.b.c. and a.p.b.c.

***********************************************************************/

#include "cl_dyn_includes.h"

/* Assumes LU preconditioned in congrad */
#ifndef LU
BOMB THE COMPILE
#endif

int pcac_x_dbl() {

   register int i,j,k,xx;
   register site *s;
   register complex cc;
   int iters=0;
   int yztvolume=nt*ny*nz;
   complex *sinkp,*sinkg1p,*sinka;         /* nx dot products         */
   Real *corrpp,*corrpa,*corrap,*corraa;   /* nx results              */
#ifdef V_TOO
   complex *sinkv,*sinkg35p;
   Real *corrvv;
#endif

   sinkp = (complex *)malloc(nx*sizeof(complex));
   sinkg1p = (complex *)malloc(nx*sizeof(complex));
   sinka = (complex *)malloc(nx*sizeof(complex));
#ifdef V_TOO
   sinkv = (complex *)malloc(nx*sizeof(complex));
   sinkg35p = (complex *)malloc(nx*sizeof(complex));
#endif
   corrpp = (Real *)malloc(nx*sizeof(Real));
   corrpa = (Real *)malloc(nx*sizeof(Real));
   corrap = (Real *)malloc(nx*sizeof(Real));
   corraa = (Real *)malloc(nx*sizeof(Real));
#ifdef V_TOO
   corrvv = (Real *)malloc(nx*sizeof(Real));
#endif

/***  generate gaussian random source on x-slice x=0 only  *************
*  normalization of random sources:        SEE libraries/gaussrand.c   *
*  <eta|eta> = <Re eta|Re eta> + <Im eta|Im eta> = 1/2 + 1/2 = 1      */

   FORALLSITES(i,s){
      for(k=0;k<4;k++)for(j=0;j<DIMF;j++){
         if(s->x == 0){
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

/* Inversions. On all but last inversion, must not pass g_rand to inverter
   because it'll be overwritten                */

/* link is initially p.b.c. in x dir, flip to a.p.b.c.      */

   boundary_flip_x(MINUS);

#ifdef V_TOO
/* compute 1/M_apbc |gamma_5*gamma_3*src>              */

   FORALLSITES(i,s) {
      mult_by_gamma( &(s->g_rand), &(s->psi[0]), ZUP );
      mult_by_gamma( &(s->psi[0]), &(s->chi[0]), GAMMAFIVE );
   }
   qic.start_flag = 0;   /* Use zero initial guess for dest */
   iters += wilson_invert(F_OFFSET(chi[0]),F_OFFSET(invv),F_OFFSET(r),
                         CURINVERT,&qic,(void *)&dcp);
#endif

/* compute 1/M_apbc |gamma_1*src>              */

   FORALLSITES(i,s) {
      mult_by_gamma( &(s->g_rand), &(s->chi[0]), XUP );
   }
   qic.start_flag = 0;   /* Use zero initial guess for dest */
   iters += wilson_invert(F_OFFSET(chi[0]),F_OFFSET(inva),F_OFFSET(r),
                         CURINVERT,&qic,(void *)&dcp);

/* compute 1/M_apbc |src>                      */

   FORALLSITES(i,s) {
      copy_wvec( &(s->g_rand), &(s->chi[0]) );
   }
   qic.start_flag = 0;   /* Use zero initial guess for dest */
   iters += wilson_invert(F_OFFSET(chi[0]),F_OFFSET(invp),F_OFFSET(r),
                         CURINVERT,&qic,(void *)&dcp);

/* link back to p.b.c.                         */
   boundary_flip_x(PLUS);

#ifdef V_TOO
/* compute 1/M_pbc |gamma_5*gamma_3*src>              */

   FORALLSITES(i,s) {
      mult_by_gamma( &(s->g_rand), &(s->psi[0]), ZUP );
      mult_by_gamma( &(s->psi[0]), &(s->chi[0]), GAMMAFIVE );
   }
   qic.start_flag = 0;   /* Use zero initial guess for dest */
   iters += wilson_invert(F_OFFSET(chi[0]),F_OFFSET(psi[0]),F_OFFSET(r),
                         CURINVERT,&qic,(void *)&dcp);

/* compute 1/M_cbc |gamma_5*gamma_3*src>              */

   FORALLSITES(i,s) {
      add_wilson_vector( &(s->psi[0]), &(s->invv), &(s->invv) );
   }
#endif

/* compute 1/M_pbc |gamma_1*src>               */

   FORALLSITES(i,s) {
      mult_by_gamma( &(s->g_rand), &(s->chi[0]), XUP );
   }
   qic.start_flag = 0;   /* Use zero initial guess for dest */
   iters += wilson_invert(F_OFFSET(chi[0]),F_OFFSET(psi[0]),F_OFFSET(r),
                         CURINVERT,&qic,(void *)&dcp);

/* compute 1/M_cbc |gamma_1*src>               */

   FORALLSITES(i,s) {
      add_wilson_vector( &(s->psi[0]), &(s->inva), &(s->inva) );
   }

/* compute 1/M_pbc |src>                       */
/* last inversion, can overwrite g_rand        */

   qic.start_flag = 0;   /* Use zero initial guess for dest */
   iters += wilson_invert(F_OFFSET(g_rand),F_OFFSET(psi[0]),F_OFFSET(r),
                         CURINVERT,&qic,(void *)&dcp);

/* compute 1/M_cbc |src>                       */

   FORALLSITES(i,s) {
      add_wilson_vector( &(s->psi[0]), &(s->invp), &(s->invp) );
   }

/***  observables  *****************************************************

all observables at py=pz=pt=0.
yztvolume = 3-volume

|src>  = random gaussian vector on x-slice = 0
|sink> = random gaussian vector on x-slice = x
|invp> = 1/M_cbc |src>
|inva> = 1/M_cbc |gamma1*src>
|invv> = 1/M_cbc |gamma_5*gamma_3*src>

A_1 ~ psibar gamma_5 gamma_1 psi,   P ~ psibar gamma_5 psi
V_3 ~ psibar gamma_3 psi

XXPP = -<pi(x) pi(0)>  = <sink|invp> <invp|sink> / yztvolume
XXAP = <A_1(x) pi(0)>  = -Re <sink*gamma_1|invp> <invp|sink> / yztvolume
XXPA = -<pi(x) A_1(0)> = -Re <sink|inva> <invp|sink> / yztvolume
XXAA = <A_1(x) A_1(0)> = Re <sink|inva> <invp|gamma_1*sink> / yztvolume
XXVV = <V_3(x) V_3(0)> = Re <sink|invv> <invp|gamma_5*gamma_3*sink> / svolume

Output format:

PROP TOTAL_ITER iters
PROP XXPP result(x=0) ... result(x=nx-1)
PROP XXAP result(x=0) ... result(x=nx-1)
PROP XXPA result(x=0) ... result(x=nx-1)
PROP XXAA result(x=0) ... result(x=nx-1)
PROP XXVV result(x=0) ... result(x=nx-1)

*/

/***  sink gaussian random vector  ************************************/

   FORALLSITES(i,s){
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
      sinkp[xx] = cmplx(0.0,0.0);
      sinkg1p[xx] = cmplx(0.0,0.0);
      sinka[xx] = cmplx(0.0,0.0);
#ifdef V_TOO
      sinkv[xx] = cmplx(0.0,0.0);
      sinkg35p[xx] = cmplx(0.0,0.0);
#endif
   }

   FORALLSITES(i,s){
      cc = wvec_dot( &(s->g_rand), &(s->invp) );
      CSUM(sinkp[s->x],cc);
      mult_by_gamma( &(s->g_rand), &(s->chi[0]), XUP );
      cc = wvec_dot( &(s->chi[0]), &(s->invp) );
      CSUM(sinkg1p[s->x],cc);
      cc = wvec_dot( &(s->g_rand), &(s->inva) );
      CSUM(sinka[s->x],cc);
#ifdef V_TOO
      mult_by_gamma( &(s->g_rand), &(s->psi[0]), ZUP );
      mult_by_gamma( &(s->psi[0]), &(s->chi[0]), GAMMAFIVE );
      cc = wvec_dot( &(s->chi[0]), &(s->invp) );
      CSUM(sinkg35p[s->x],cc);
      cc = wvec_dot( &(s->g_rand), &(s->invv) );
      CSUM(sinkv[s->x],cc);
#endif
   }

   for(xx=0;xx<nx;xx++){
      g_complexsum( &(sinkp[xx]) );
      g_complexsum( &(sinkg1p[xx]) );
      g_complexsum( &(sinka[xx]) );
#ifdef V_TOO
      g_complexsum( &(sinkv[xx]) );
      g_complexsum( &(sinkg35p[xx]) );
#endif
   }

/***  compute correlators  ********************************************/

   for(xx=0;xx<nx;xx++){
      corrpp[xx] = sinkp[xx].real*sinkp[xx].real
                   +sinkp[xx].imag*sinkp[xx].imag;
      corrap[xx] = sinkg1p[xx].real*sinkp[xx].real
                   +sinkg1p[xx].imag*sinkp[xx].imag;
      corrpa[xx] = sinka[xx].real*sinkp[xx].real
                   +sinka[xx].imag*sinkp[xx].imag;
      corraa[xx] = sinka[xx].real*sinkg1p[xx].real
                   +sinka[xx].imag*sinkg1p[xx].imag;
#ifdef V_TOO
      corrvv[xx] = sinkv[xx].real*sinkg35p[xx].real
                   +sinkv[xx].imag*sinkg35p[xx].imag;
#endif

/* normalization and sign                     */

      corrpp[xx] *= (4*kappa*kappa)/yztvolume;
      corrap[xx] *= (-4*kappa*kappa)/yztvolume;
      corrpa[xx] *= (-4*kappa*kappa)/yztvolume;
      corraa[xx] *= (4*kappa*kappa)/yztvolume;
#ifdef V_TOO
      corrvv[xx] *= (-4*kappa*kappa)/yztvolume;
#endif
   }

/***  printouts  ******************************************************/

   if(this_node==0) {
      printf("PROP TOTAL_ITERS %d\n", iters);

/* print x-slices correlations                */

      printf("PROP XXPP");
      for(xx=0;xx<nx;xx++){
         printf(" %e", (double)corrpp[xx]);
      }
      printf("\n");

      printf("PROP XXAP");
      for(xx=0;xx<nx;xx++){
         printf(" %e", (double)corrap[xx]);
      }
      printf("\n");

      printf("PROP XXPA");
      for(xx=0;xx<nx;xx++){
         printf(" %e", (double)corrpa[xx]);
      }
      printf("\n");

      printf("PROP XXAA");
      for(xx=0;xx<nx;xx++){
         printf(" %e", (double)corraa[xx]);
      }
      printf("\n");

#ifdef V_TOO
      printf("PROP XXVV");
      for(xx=0;xx<nx;xx++){
         printf(" %e", (double)corrvv[xx]);
      }
      printf("\n");
#endif

   }

   free(sinkp);  free(sinkg1p);  free(sinka);
   free(corrpp);  free(corrpa);  free(corrap);  free(corraa);
#ifdef V_TOO
     free(sinkv);  free(sinkg35p);  free(corrvv);
#endif

   fflush(stdout);
   return(iters);
}
/* end pcac_x_dbl */

