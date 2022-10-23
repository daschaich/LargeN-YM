/* Fix Coulomb or Lorentz gauge by doing successive SU(2) gauge hits */
/* Uses double precision global sums */
/* This version does automatic reunitarization at preset intervals */

/* This version for arbitrary NCOL and does gauge fixing via relaxation in full SU(N) group */

/* development version does not use field_offset diffmat, field_offset sumvec */
/* Prototype...

   void gaugefix(int gauge_dir,Real relax_boost,int max_gauge_iter,
        Real gauge_fix_tol, field_offset diffmat, field_offset sumvec,
        int nvector, field_offset vector_offset[], int vector_parity[],
        int nantiherm, field_offset antiherm_offset[],
        int antiherm_parity[] )
   -------------------------------------------------------------------

   NOTE: For staggered fermion applications, it is necessary to remove
   the KS phases from the gauge links before calling this procedure.
   See "rephase" in setup.c.

   -------------------------------------------------------------------
   EXAMPLE:  Fixing only the link matrices to Coulomb gauge with scratch
     space in mp (matrix) and chi (vector):

   gaugefix(TUP,(Real)1.5,500,(Real)1.0e-7,
         F_OFFSET(mp),F_OFFSET(chi),0,NULL,NULL,0,NULL,NULL);

   -------------------------------------------------------------------
   EXAMPLE:  Fixing Coulomb gauge with respect to the y direction
      in the staggered fermion scheme and simultaneously transforming
      the pseudofermion fields and gauge-momenta involved in updating:

   int nvector = 3;
   field_offset vector_offset[3] = { F_OFFSET(g_rand), F_OFFSET(phi),
        F_OFFSET(xxx) };
   int vector_parity[3] = { EVENANDODD, EVEN, EVEN };
   int nantiherm = 4;
   field_offset antiherm_offset[4] = { F_OFFSET(mom[0]), F_OFFSET(mom[1]),
       F_OFFSET(mom[2]), F_OFFSET(mom[3]) };
   field_offset antiherm_parity[4] = { EVENANDODD, EVENANDODD, EVENANDODD,
       EVENANDODD }

   rephase( OFF );
   gaugefix(YUP,(Real)1.8,500,(Real)2.0e-6,
       F_OFFSET(tempmat1),F_OFFSET(tempvec[0]),
       nvector,vector_offset,vector_parity,
       nantiherm,antiherm_offset,antiherm_parity);
   rephase( ON );

   -------------------------------------------------------------------

   gauge_dir     specifies the direction of the "time"-like hyperplane
                 for the purposes of defining Coulomb or Lorentz gauge
      TUP    for evaluating propagators in the time-like direction
      ZUP    for screening lengths.
      8      for Lorentz gauge
   relax_boost     Overrelaxation parameter
   max_gauge_iter  Maximum number of iterations
   gauge_fix_tol   Stop if change is less than this
   diffmat         Scratch space for an su3 matrix
   sumvec          Scratch space for an su3 vector
   NOTE: if diffmat or sumvec are negative, gaugefix mallocs its own
   scratch space. */

#include "generic_includes.h"
#include "../include/generic_nhyp.h"
#define REUNIT_INTERVAL 10
#define EPS_SQ 1.0e-8

/* Scratch space */
static matrix *diffmatp;               /* malloced diffmat pointer */
static vector *sumvecp;                /* malloced sumvec pointer */
field_offset diffmat_offset,sumvec_offset;  /* field offsets */

void gaugefix(int gauge_dir,Real relax_boost,int max_gauge_iter,
              Real gauge_fix_tol, field_offset diffmat, field_offset sumvec,
              int nvector, field_offset vector_offset[], int vector_parity[],
              int nantiherm, field_offset antiherm_offset[],
              int antiherm_parity[]) {

  int gauge_iter;
  double current_av, old_av = 0, del_av = 0;
  void gaugefixscratch(field_offset diffmat, field_offset sumvec);
  void gaugefixstep(int gauge_dir,double *av_gauge_fix_action,Real relax_boost,
        int nvector, field_offset vector_offset[], int vector_parity[],
        int nantiherm, field_offset antiherm_offset[],
        int antiherm_parity[]);

  // We require at least 8 gen_pt values for gauge fixing
  if(N_POINTERS < 8) {
      printf("gaugefix: N_POINTERS must be at least %d.  Fix the code.\n",
       N_POINTERS);
      fflush(stdout); terminate(1);
    }


  /* Set up work space */
  gaugefixscratch(diffmat,sumvec);

  /* Do at most max_gauge_iter iterations, but stop after the second step if */
  /* the change in the avg gauge fixing action is smaller than gauge_fix_tol */

  for (gauge_iter=0; gauge_iter < max_gauge_iter; gauge_iter++)
    {
      gaugefixstep(gauge_dir,&current_av,relax_boost,
       nvector, vector_offset, vector_parity,
       nantiherm, antiherm_offset, antiherm_parity);

      if(gauge_iter != 0)
  {
    del_av = current_av - old_av;
    if (fabs(del_av) < gauge_fix_tol) break;
  }
      old_av = current_av;

      /* Reunitarize when iteration count is a multiple of REUNIT_INTERVAL */
      if((gauge_iter % REUNIT_INTERVAL) == (REUNIT_INTERVAL - 1))
  {
    node0_printf("step %d av gf action %.8e, delta %.3e\n",
           gauge_iter,current_av,del_av);
    reunitarize();
  }
    }
  /* Reunitarize at the end, unless we just did it in the loop */
  if((gauge_iter % REUNIT_INTERVAL) != 0)
    reunitarize();

  /* Free workspace */
  free(diffmatp);
  free(sumvecp);

  if(this_node==0)
    printf("GFIX: Ended at step %d. Av gf action %.8e, delta %.3e\n",
     gauge_iter,(double)current_av,(double)del_av);
}

void gaugefixstep(int gauge_dir,double *av_gauge_fix_action,Real relax_boost,
      int nvector, field_offset vector_offset[], int vector_parity[],
      int nantiherm, field_offset antiherm_offset[],
      int antiherm_parity[] )
{
  /* Carry out one iteration in the gauge-fixing process */

  double get_gauge_fix_action(int gauge_dir,int parity);
  void do_hit_full(int gauge_dir, int parity,  Real relax_boost,
       int nvector, field_offset vector_offset[], int vector_parity[],
       int nantiherm, field_offset antiherm_offset[],
       int antiherm_parity[] );
  int parity;
  msg_tag *mtag[8];
  Real gauge_fix_action;
  register int dir,i;
  register site *s;

  /* Alternate parity to prevent interactions during gauge transformation */
  *av_gauge_fix_action = 0.;
  g_sync();
  fflush(stdout);

  for(parity = ODD; parity <= EVEN; parity++)
    {
      /* Start gathers of downward links */

      FORALLUPDIR(dir)
  {
    mtag[dir] = start_gather_site( F_OFFSET(linkf[dir]), sizeof(matrix),
           OPP_DIR(dir), parity, gen_pt[dir] );
  }

      /* Wait for gathers */

      FORALLUPDIR(dir)
  {
    wait_gather(mtag[dir]);
  }

      /* Total gauge fixing action for sites of this parity: Before */
      gauge_fix_action = get_gauge_fix_action(gauge_dir,parity);

      /* Do optimum gauge hit on various subspaces */

      do_hit_full(gauge_dir,parity, relax_boost,
      nvector, vector_offset, vector_parity,
      nantiherm, antiherm_offset, antiherm_parity);
      /* Total gauge fixing action for sites of this parity: After */
      gauge_fix_action = get_gauge_fix_action(gauge_dir,parity);

      *av_gauge_fix_action += gauge_fix_action;

      /* Scatter downward link matrices by gathering to sites of */
      /* opposite parity */

      FORALLUPDIR(dir)
  {
    /* Synchronize before scattering to be sure the new modified link */
    /* matrices are all ready to be scattered and diffmat is not */
    /* overwritten before it is used */
    g_sync();

    /* First copy modified link for this dir */
    /* from comm buffer or node to diffmat */

    FORSOMEPARITY(i,s,parity)
      {
        mat_copy_f((matrix *)(gen_pt[dir][i]), &diffmatp[i]);
      }

    /* Now we are finished with gen_pt[dir] */
    cleanup_gather(mtag[dir]);

    /* Synchronize to make sure the previous copy happens before the */
    /* subsequent gather below  */
    g_sync();

    /* Gather diffmat onto sites of opposite parity */

    mtag[dir] = start_gather_field( diffmatp, sizeof(matrix),
            dir, OPP_PAR(parity), gen_pt[dir] );

    wait_gather(mtag[dir]);

    /* Copy modified matrices into proper location */

    FORSOMEPARITY(i,s,OPP_PAR(parity))
      mat_copy_f((matrix *)(gen_pt[dir][i]),&(s->linkf[dir]));

    cleanup_gather(mtag[dir]);
  }

    }
} /* gaugefixstep */

void do_hit_full(int gauge_dir, int parity, Real relax_boost,
     int nvector, field_offset vector_offset[], int vector_parity[],
     int nantiherm, field_offset antiherm_offset[],
     int antiherm_parity[] )
{
  /* Do relaxation in SU(NCOL) */
  register int dir,i;
  register site *s;
  matrix u;
  void accum_gauge_hit(int gauge_dir,int parity);
  complex find_det(matrix *c);
  Real f[NCOL],theta;
  complex tt1,phase,det;
  int j,k;
  matrix  Omega, Q[NCOL];
#if (NCOL>2)
  matrix eQ;
#endif


  /* Accumulate sums for determining optimum gauge hit */

  accum_gauge_hit(gauge_dir,parity);

  FORSOMEPARITYDOMAIN(i,s,parity)
    {
    /* begin general code, stolen from nhyp */
      mat_copy_f(&diffmatp[i],&Omega);
      mult_an_f(&Omega,&Omega,&Q[1]);
#ifndef NHYP_DEBUG
      compute_fhb(&Q[1],f,NULL, 0);
#else
      compute_fhb(&Omega,&Q[1],f,NULL, 0);
#endif
#if (NCOL==2)
      scalar_mult_mat_f(&Omega,f[0],&u);
#else
      /* compute Q^(-1/2) via Eq. (3.8)  */
      scalar_mult_mat_f(&Q[1],f[1],&eQ);

      for(j=2;j<NCOL;j++){
  mult_nn_f(&Q[1],&Q[j-1],&Q[j]);
  scalar_mult_add_mat_f(&eQ,&Q[j],f[j],&eQ);
      }
      scalar_add_diag_f(&eQ,f[0]);

      /* multiply Omega by eQ = (Omega^\dagger Omega)^(-1/2) to make u in U(N) */
      mult_nn_f(&Omega,&eQ,&u);
#endif

      /* check det...*/
      det=find_det(&u);

      /* spread the phase over everybody */
      theta=atan(det.imag/(det.real+EPS_SQ));
      phase=ce_itheta(-theta/NCOL);

      for(j=0;j<NCOL;j++)for(k=0;k<NCOL;k++){
  CMUL(phase,u.e[j][k],tt1);
  u.e[j][k]=tt1;
      }

      /* Do gauge transformation on all upward links */

      FORALLUPDIR(dir){
  mult_nn_f(&u,&(s->linkf[dir]),&Omega);
  mat_copy_f(&Omega,&(s->linkf[dir]));
      }

      /* Do gauge transformation hit on all downward links */

      FORALLUPDIR(dir){
  mult_na_f((matrix *)gen_pt[dir][i],&u,&Omega);
  mat_copy_f(&Omega,(matrix *)gen_pt[dir][i]);
      }
      /* code to transform other variables has been removed */
    }
  /* Exit with modified downward links left in communications buffer */
}
/* do_hit */

void accum_gauge_hit(int gauge_dir,int parity)
{

  /* Accumulates environment of link matrices for  gauge fixing, result in diffmatp */


  /* To be general, we want to maximize
     Re Tr V(x) \Sigma^\dagger where \Sigma^\datter \sum_nu(U_\nu(x)+U_\nu(x-\nu)^\dagger
     and our recycled nHYP code finds the U(N) matrix V ~ Sigma*/

  register int j,k;
  register matrix *m1,*m2;
  register int dir,i;
  register site *s;
  matrix m11;

  /* Clear diffmat */

  FORSOMEPARITY(i,s,parity)
  clear_mat_f(&diffmatp[i]);

  /* Subtract upward link contributions */

  FORSOMEPARITYDOMAIN(i,s,parity) {
      FORALLUPDIRBUT(gauge_dir,dir) {
    /* Upward link matrix */
    m1 = &(s->linkf[dir]);
    for(j=0;j<NCOL;j++)for(k=0;k<NCOL;k++){
      m11.e[j][k]  = conjg(&(m1->e[k][j]));
    }
    add_mat_f( &diffmatp[i], &m11, &diffmatp[i]);
  }
    }

  /* Add downward link contributions */
  FORSOMEPARITYDOMAIN(i,s,parity) {
      FORALLUPDIRBUT(gauge_dir,dir) {
    /* Downward link matrix */
    m2 = (matrix *)gen_pt[dir][i];
    add_mat_f( &diffmatp[i], m2, &diffmatp[i]);
  }
    }
} /* accum_gauge_hit */




double get_gauge_fix_action(int gauge_dir,int parity)
{
  /* Adds up the gauge fixing action for sites of given parity */
  /* Returns average over these sites */
  /* The average is normalized to a maximum of 1 when all */
  /* links are unit matrices */

  register int dir,i,ndir;
  register site *s;
  register matrix *m1, *m2;
  double gauge_fix_action;
  complex trace;

  gauge_fix_action = 0.0;

  FORSOMEPARITYDOMAIN(i,s,parity)
    {
      FORALLUPDIRBUT(gauge_dir,dir)
  {
    m1 = &(s->linkf[dir]);
    m2 = (matrix *)gen_pt[dir][i];

    trace = trace_f(m1);
    gauge_fix_action += (double)trace.real;

    trace = trace_f(m2);
    gauge_fix_action += (double)trace.real;
  }
    }

  /* Count number of terms to average */
  ndir = 0; FORALLUPDIRBUT(gauge_dir,dir)ndir++;

  /* Sum over all sites of this parity */
  g_doublesum( &gauge_fix_action);

  /* Average is normalized to max of 1/2 on sites of one parity */
  return(gauge_fix_action /((double)(2*NCOL*ndir*nx*ny*nz*nt)));
} /* get_gauge_fix_action */

void gaugefixscratch(field_offset diffmat, field_offset sumvec) {
  diffmat_offset = diffmat;
  diffmatp = NULL;

  diffmatp = (matrix *)malloc(sizeof(matrix)*sites_on_node);
  if(diffmatp == NULL)
    {
      node0_printf("gaugefix: Can't malloc diffmat\n");
      fflush(stdout);terminate(1);
    }
} /* gaugefixscratch */

