/*
* wilson/clover fermions

* Fermionic observables to find kappa_c by imposing PCAC relation:

* 4 = number of Dirac components
* NCOL = number of colors
* DIMF = size of fermion rep

* clover operator:
* M = A - kappa*( Dslash_eo + DSLASH_oe )
* 2*kappa=1/(4+m0)

* measurements (vec{p}=0):
* TTPP = -<pi(t) pi(0)>
* TTAP = <A_4(t) pi(0)>
* TTPA = -<pi(t) A_4(0)>
* TTAA = <A_4(t) A_4(0)>

if V_TOO is set in defines.h we also do:
* TTVV = <V_3(t) V_3(0)>

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

int pcac_t() {

   register int i,j,k,tt;
   register site *s;
   register complex cc;
   int iters=0;
   int svolume=nx*ny*nz;
   complex *sinkp,*sinkg4p,*sinka;         /* nt dot products  */
   Real *corrpp,*corrpa,*corrap,*corraa;   /* nt results       */
#ifdef V_TOO
   complex *sinkv,*sinkg35p;
   Real *corrvv;
#endif

   sinkp = (complex *)malloc(nt*sizeof(complex));
   sinkg4p = (complex *)malloc(nt*sizeof(complex));
   sinka = (complex *)malloc(nt*sizeof(complex));
#ifdef V_TOO
   sinkv = (complex *)malloc(nt*sizeof(complex));
   sinkg35p = (complex *)malloc(nt*sizeof(complex));
#endif
   corrpp = (Real *)malloc(nt*sizeof(Real));
   corrpa = (Real *)malloc(nt*sizeof(Real));
   corrap = (Real *)malloc(nt*sizeof(Real));
   corraa = (Real *)malloc(nt*sizeof(Real));
#ifdef V_TOO
   corrvv = (Real *)malloc(nt*sizeof(Real));
#endif

/***  generate gaussian random source on time-slice t=0 only  **********
*  normalization of random sources:        SEE libraries/gaussrand.c   *
*  <eta|eta> = <Re eta|Re eta> + <Im eta|Im eta> = 1/2 + 1/2 = 1      */

   FORALLSITES(i,s){
      for(k=0;k<4;k++)for(j=0;j<DIMF;j++){
         if(s->t == 0){
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

#ifdef V_TOO
/* compute 1/M_apbc |gamma_5*gamma_3*src>              */

   FORALLSITES(i,s) {
      mult_by_gamma( &(s->g_rand), &(s->psi[0]), ZUP );
      mult_by_gamma( &(s->psi[0]), &(s->chi[0]), GAMMAFIVE );
   }
   qic.start_flag = 0;   /* Use zero initial guess for dest */
   iters += wilson_invert(F_OFFSET(chi[0]),F_OFFSET(psi[1]),F_OFFSET(r),
                         CURINVERT,&qic,(void *)&dcp);
#endif

/* compute 1/M_apbc |gamma_4*src>              */

   FORALLSITES(i,s) {
      mult_by_gamma( &(s->g_rand), &(s->chi[0]), TUP );
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

/* link is initially a.p.b.c., flip to p.b.c.  */
   boundary_flip(PLUS);

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
      add_wilson_vector( &(s->psi[0]), &(s->psi[1]), &(s->psi[1]) );
   }
#endif

/* compute 1/M_pbc |gamma_4*src>               */

   FORALLSITES(i,s) {
      mult_by_gamma( &(s->g_rand), &(s->chi[0]), TUP );
   }
   qic.start_flag = 0;   /* Use zero initial guess for dest */
   iters += wilson_invert(F_OFFSET(chi[0]),F_OFFSET(psi[0]),F_OFFSET(r),
                         CURINVERT,&qic,(void *)&dcp);

/* compute 1/M_cbc |gamma_4*src>               */

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

/* link back to a.p.b.c.                       */
   boundary_flip(MINUS);

/***  observables  *****************************************************

all observables at vec{p}=0.
svolume = 3-volume

|src>  = random gaussian vector on time-slice = 0
|sink> = random gaussian vector on time-slice = t
|invp> = 1/M_cbc |src>
|inva> = 1/M_cbc |gamma_4*src>
|invv> = 1/M_cbc |gamma_5*gamma_3*src>
NOTE: to save space, we use psi[1] for invv.

A_4 ~ psibar gamma_5 gamma_4 psi,   P ~ psibar gamma_5 psi
V_3 ~ psibar gamma_3 psi

TTPP = -<pi(t) pi(0)>  = <sink|invp> <invp|sink> / svolume
TTAP = <A_4(t) pi(0)>  = -Re <sink*gamma_4|invp> <invp|sink> / svolume
TTPA = -<pi(t) A_4(0)> = -Re <sink|inva> <invp|sink> / svolume
TTAA = <A_4(t) A_4(0)> = Re <sink|inva> <invp|gamma_4*sink> / svolume
TTVV = <V_3(t) V_3(0)> = Re <sink|invv> <invp|gamma_5*gamma_3*sink> / svolume

Output format:

PROP TOTAL_ITER iters
PROP TTPP result(t=0) ... result(t=nt-1)
PROP TTAP result(t=0) ... result(t=nt-1)
PROP TTPA result(t=0) ... result(t=nt-1)
PROP TTAA result(t=0) ... result(t=nt-1)
PROP TTVV result(t=0) ... result(t=nt-1)

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

   for(tt=0;tt<nt;tt++){
      sinkp[tt] = cmplx(0.0,0.0);
      sinkg4p[tt] = cmplx(0.0,0.0);
      sinka[tt] = cmplx(0.0,0.0);
#ifdef V_TOO
      sinkv[tt] = cmplx(0.0,0.0);
      sinkg35p[tt] = cmplx(0.0,0.0);
#endif
   }

   FORALLSITES(i,s){
      cc = wvec_dot( &(s->g_rand), &(s->invp) );
      CSUM(sinkp[s->t],cc);
      mult_by_gamma( &(s->g_rand), &(s->chi[0]), TUP );
      cc = wvec_dot( &(s->chi[0]), &(s->invp) );
      CSUM(sinkg4p[s->t],cc);
      cc = wvec_dot( &(s->g_rand), &(s->inva) );
      CSUM(sinka[s->t],cc);
#ifdef V_TOO
      mult_by_gamma( &(s->g_rand), &(s->psi[0]), ZUP );
      mult_by_gamma( &(s->psi[0]), &(s->chi[0]), GAMMAFIVE );
      cc = wvec_dot( &(s->chi[0]), &(s->invp) );
      CSUM(sinkg35p[s->t],cc);
      cc = wvec_dot( &(s->g_rand), &(s->psi[1]) );
      CSUM(sinkv[s->t],cc);
#endif
   }

   for(tt=0;tt<nt;tt++){
      g_complexsum( &(sinkp[tt]) );
      g_complexsum( &(sinkg4p[tt]) );
      g_complexsum( &(sinka[tt]) );
#ifdef V_TOO
      g_complexsum( &(sinkv[tt]) );
      g_complexsum( &(sinkg35p[tt]) );
#endif
   }

/***  compute correlators  ********************************************/

   for(tt=0;tt<nt;tt++){
      corrpp[tt] = sinkp[tt].real*sinkp[tt].real
                   +sinkp[tt].imag*sinkp[tt].imag;
      corrap[tt] = sinkg4p[tt].real*sinkp[tt].real
                   +sinkg4p[tt].imag*sinkp[tt].imag;
      corrpa[tt] = sinka[tt].real*sinkp[tt].real
                   +sinka[tt].imag*sinkp[tt].imag;
      corraa[tt] = sinka[tt].real*sinkg4p[tt].real
                   +sinka[tt].imag*sinkg4p[tt].imag;
#ifdef V_TOO
      corrvv[tt] = sinkv[tt].real*sinkg35p[tt].real
                   +sinkv[tt].imag*sinkg35p[tt].imag;
#endif

/* normalization and sign                     */

      corrpp[tt] *= (4*kappa*kappa)/svolume;
      corrap[tt] *= (-4*kappa*kappa)/svolume;
      corrpa[tt] *= (-4*kappa*kappa)/svolume;
      corraa[tt] *= (4*kappa*kappa)/svolume;
#ifdef V_TOO
      corrvv[tt] *= (-4*kappa*kappa)/svolume;
#endif
   }

/***  printouts  ******************************************************/

   if(this_node==0) {
      printf("PROP TOTAL_ITERS %d\n", iters);

/* print time-slices correlations             */

      printf("PROP TTPP");
      for(tt=0;tt<nt;tt++){
         printf(" %e", (double)corrpp[tt]);
      }
      printf("\n");

      printf("PROP TTAP");
      for(tt=0;tt<nt;tt++){
         printf(" %e", (double)corrap[tt]);
      }
      printf("\n");

      printf("PROP TTPA");
      for(tt=0;tt<nt;tt++){
         printf(" %e", (double)corrpa[tt]);
      }
      printf("\n");

      printf("PROP TTAA");
      for(tt=0;tt<nt;tt++){
         printf(" %e", (double)corraa[tt]);
      }
      printf("\n");

#ifdef V_TOO
      printf("PROP TTVV");
      for(tt=0;tt<nt;tt++){
         printf(" %e", (double)corrvv[tt]);
      }
      printf("\n");
#endif

   }

   free(sinkp);  free(sinkg4p);  free(sinka);
   free(corrpp);  free(corrpa);  free(corrap);  free(corraa);
#ifdef V_TOO
     free(sinkv);  free(sinkg35p);  free(corrvv);
#endif

   fflush(stdout);
   return(iters);
}
/* end pcac_t */

