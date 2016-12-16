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
int meascount,todo;
int prompt;
double dssplaq,dstplaq,dssplaq_frep,dstplaq_frep;
int m_iters=0,s_iters,avm_iters,avs_iters;
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
#ifdef TIMING
        time_dcongrad=time_fermion_force=time_fermion_rep=time_block_nhyp=time_compute_fhb=0.;
// time_gmp=time_jacobi=0.;
#endif

	/* call plaquette measuring process */
/* Check: compute initial plaquette (T. D.) */
                d_plaquette(&dssplaq,&dstplaq);
                if(this_node==0)printf("START %e %e %e\n",
                    dssplaq,dstplaq,
                        dssplaq+dstplaq);

	for(todo=warms; todo > 0; --todo ){
	    update();
	}
	if(this_node==0)printf("WARMUPS COMPLETED\n");

	/* perform measuring trajectories, reunitarizing and measuring 	*/
	meascount=0;		/* number of measurements 		*/
	plp = cmplx((Real)99.9,(Real)99.9);
	avm_iters = avs_iters = 0;
#ifdef SPECTRUM
        avspect_iters = 0;
#endif
	for(todo=trajecs; todo > 0; --todo ){
	    /* do the trajectories */
	    s_iters=update();

	    /* Do "local" measurements every trajectory! */
#ifdef SF   /* compute the SF coupling */
            coupling();
#endif
            /* The action from the RG trans */
#ifdef IMP  /* For improved action only.
               Plaqutte action trivially follows from ds[s|t]plaq */
            gauge_action(&dssplaq);
            if(this_node==0)printf("ACTION_V  %e  %e\n",
                dssplaq,(dssplaq)/(double)(volume*6));
#endif
	    /* call the Polyakov loop measuring program */
	    plp = ploop();

	    /* call plaquette measuring process */
	    d_plaquette(&dssplaq,&dstplaq);
	    d_plaquette_frep(&dssplaq_frep,&dstplaq_frep);
#ifdef LU
            /* generate a pseudofermion configuration */
	    m_iters = f_measure_cl();
#endif
	    ++meascount;
	    avm_iters += m_iters;
	    avs_iters += s_iters;

	    if(this_node==0)printf("GMES %e %e %e %e %e %e %e\n",
		(double)plp.real,(double)plp.imag,(double)m_iters,
		dssplaq,dstplaq,dssplaq_frep,dstplaq_frep);
	    /* Re(Polyakov) Im(Poyakov) cg_iters ss_plaq st_plaq */

	    /* measure other stuff every "propinterval" trajectories */
	    if(((todo-1)%propinterval) == 0){
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

	      fixflag = NO_GAUGE_FIX;
#ifdef SPECTRUM
#ifdef SCREEN
		gaugefix(ZUP,(Real)1.5,500,(Real)GAUGE_FIX_TOL,
			 F_OFFSET(staple),F_OFFSET(tempmat1),
			 0,NULL,NULL,0,NULL,NULL);
		fixflag = COULOMB_GAUGE_FIX;
		spect_iters = s_props_cl();
		avspect_iters += spect_iters;
#else	/* spectrum in time direction */
		gaugefix(TUP,(Real)1.5,500,(Real)GAUGE_FIX_TOL,
			 F_OFFSET(staple),F_OFFSET(tempmat1),
			 0,NULL,NULL,0,NULL,NULL);
		fixflag = COULOMB_GAUGE_FIX;
/* commented 15 OCT 95, MBW....we'll use periodic b.c. for spect. */
/* Don't need t_props if we are doing w_spectrum  - C. DeTar */
/*		spect_iters = t_props_cl();
		avspect_iters += spect_iters; */
		spect_iters = w_spectrum_cl();
		avspect_iters += spect_iters;
#endif	/* end ifndef SCREEN */
#endif	/* end ifdef SPECTRUM */

	    }
	    fflush(stdout);

	}	/* end loop over trajectories */

	if(this_node==0)printf("RUNNING COMPLETED\n");
/* Check: compute final plaquette (T. D.) */
                d_plaquette(&dssplaq,&dstplaq);
                if(this_node==0)printf("STOP %e %e %e\n",
                    dssplaq,dstplaq,
                        dssplaq+dstplaq);

	if(meascount>0)  {
	    if(this_node==0)printf("average cg iters for update= %e\n",
		(double)avs_iters/meascount);
	    if(this_node==0)printf("average cg iters for measurement= %e\n",
		(double)avm_iters/meascount);
#ifdef SPECTRUM
	    if(this_node==0)printf("average cg iters for spectrum = %e\n",
		(double)avspect_iters/meascount);
#endif
	}

	dtime += dclock();
	if(this_node==0){
	    printf("\nTotal time         = %e seconds\n",dtime);
#ifdef TIMING
            printf("dcongrad time      = %e seconds\n",time_dcongrad);
            printf("fermion_force time = %e seconds\n",time_fermion_force);
//          printf("fermion_rep time   = %e seconds\n",time_fermion_rep);
            printf("block_nhyp time    = %e seconds\n",time_block_nhyp);
            printf("compute_fhb time   = %e seconds\n",time_compute_fhb);
//          printf("gmp time           = %e seconds\n",time_gmp);
//          printf("jacobi time        = %e seconds\n\n",time_jacobi);
            printf("\n");
#endif
	    printf("total_iters = %d\n",total_iters);
	}
	fflush(stdout);
	dtime = -dclock();

	/* save lattice if requested */
        if( saveflag != FORGET ){
	  save_lattice( saveflag, savefile, stringLFN );
        }
    }
    return 0;
}
