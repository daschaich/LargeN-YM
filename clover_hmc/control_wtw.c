/******** control.c ********/
/* Main procedure for SU3 with dynamical clover fermions
   and Symanzik/Tadpole improved gauge action	*/
/* MIMD version 6 */

/* Modifications
   1996 Created by Matt Wingate
   6/6/98 Version 5
*/

/* This version combines code for the PHI algorithm (approriate for 4
   flavors) and the R algorithm for "epsilon squared" updating of
   1 to 4 flavors.  Compilation should occur with PHI_ALGORITHM defined
   for the former and not defined for the latter.  It also contains code
   for the hybrid Monte Carlo algorithm, for which HMC_ALGORITHM and
   PHI_ALGORITHM should be defined.  (Actually, the
   changes to control.c are minimal and the real differences will appear
   in update.c */

#define CONTROL

#include "cl_dyn_includes.h"

int main(int argc, char *argv[])  {
int meascount;
int prompt;
int m_iters=0,avm_iters,avs_iters;
#ifdef SPECTRUM
int spect_iters,avspect_iters;
#endif
complex plp;
double dtime;

setlinebuf(stdout); /*DEBUG*/
 initialize_machine(&argc,&argv);
#ifdef HAVE_QDP
  QDP_initialize(&argc, &argv);
#endif
  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);

 g_sync();
    /* set up */
    prompt = setup();

    /* loop over input sets */
    while( readin(prompt) == 0){

/* set up loop tables */
#ifdef IMP
        make_loop_table2();
#endif
	/* perform warmup trajectories */
	dtime = -dclock();


	/* perform measuring trajectories, reunitarizing and measuring 	*/
	meascount=0;		/* number of measurements 		*/
	plp = cmplx((Real)99.9,(Real)99.9);
	avm_iters = avs_iters = 0;
#ifdef SPECTRUM
        avspect_iters = 0;
#endif
	    /* measure other stuff every "propinterval" trajectories */
#ifdef LU
#ifdef PCAC
#ifndef SF
            /* correlators in pcac relation */
/*  t direction: doubling the lattice via PBC+APBC  */
            m_iters += pcac_t();
/*  x direction:  PBC and PBC+APBC                  */
            if(nt<nx) m_iters += pcac_x();
/*  x direction: PBC+APBC only
            if(nt<nx) m_iters += pcac_x_dbl();      */
#else  /* ifdef SF */
            if(nt>=4){
/*  t direction, SFBC                               */
              m_iters += pcac_sf();
/*  x direction:  PBC and PBC+APBC                  */
              if(nt<nx) m_iters += pcac_x();
            }
#endif /* SF   */
#endif /* PCAC */
#endif /* LU   */
	}
	    fflush(stdout);


	if(this_node==0)printf("RUNNING COMPLETED\n");

	dtime += dclock();
	if(this_node==0){
	    printf("Time = %e seconds\n",dtime);
	    printf("total_iters = %d\n",total_iters);
	}
	fflush(stdout);
	dtime = -dclock();

	/* save lattice if requested */
        if( saveflag != FORGET ){
	  save_lattice( saveflag, savefile, stringLFN );
        }
    
    return 0;
}
