// -----------------------------------------------------------------
// Main procedure for pure-gauge evolution
// Removed hybrid Monte Carlo updates,
// keeping just over-relaxed quasi-heat bath
#define CONTROL
#include "pg_includes.h"

int main(int argc, char *argv[]) {
  int meascount,todo;
  int prompt;
  double dssplaq,dstplaq;
  complex plp;
  double dtime;

  initialize_machine(&argc,&argv);

  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)
    terminate(1);

  g_sync();
  /* set up */
  prompt = setup();

  /* loop over input sets */
  while( readin(prompt) == 0) {

    /* perform warmup trajectories */
    dtime = -dclock();

    for(todo=warms; todo > 0; --todo)
      update();

    if(this_node==0)printf("WARMUPS COMPLETED\n");

    /* perform measuring trajectories, reunitarizing and measuring  */
    meascount=0;            /* number of measurements               */
    plp = cmplx(99.9, 99.9);
    for(todo=trajecs; todo > 0; --todo ) {

      /* do the trajectories */
      update();

      /* measure every "propinterval" trajectories */
      if((todo%propinterval) == 0){
        /* call plaquette measuring process */
        d_plaquette(&dssplaq,&dstplaq);

        /* call the Polyakov loop measuring program */
        plp = ploop();
        ++meascount;
        if(this_node==0)printf("GMES %e %e %e %e %e\n",
            (double)plp.real,(double)plp.imag,99.9,dssplaq,dstplaq);
        /* Re(Polyakov) Im(Poyakov) cg_iters ss_plaq st_plaq */
        fflush(stdout);
      }
    }       /* end loop over trajectories */

#ifdef ORA_ALGORITHM
    /* gaugefix if requested */
    if( fixflag == COULOMB_GAUGE_FIX){
      gaugefix(TUP,(Real)1.8,600,(Real)GAUGE_FIX_TOL);
      if(this_node==0)printf("FIXED TO COULOMB GAUGE\n");
      fflush(stdout);
    }
    else if( fixflag == LANDAU_GAUGE_FIX){
      gaugefix(8,(Real)1.8,600,(Real)GAUGE_FIX_TOL);
      if(this_node==0)printf("FIXED TO LANDAU GAUGE\n");
      fflush(stdout);
    }
#endif

    if(this_node==0)printf("RUNNING COMPLETED\n");

    dtime += dclock();
    if(this_node==0){
      printf("Time = %e seconds\n",dtime);
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
// -----------------------------------------------------------------
