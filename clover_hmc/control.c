// -----------------------------------------------------------------
// Main procedure for nHYP-smeared Wilson-clover SU(N) evolation
#define CONTROL
#include "cl_dyn_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int main(int argc, char *argv[]) {
  int traj_done, m_iters, s_iters, avm_iters = 0, avs_iters = 0, Nmeas = 0;
  int prompt;
  Real f_eps0, f_eps1, g_eps;
  double ss_plaq, st_plaq, ss_plaq_frep, st_plaq_frep, dtime;
#ifdef SPECTRUM
  int spect_iters, avspect_iters = 0;
#endif
  complex plp = cmplx(99.0, 99.0);

  // Set up
  setlinebuf(stdout); // DEBUG
  // Remap standard I/O
  if (remap_stdio_from_args(argc, argv) == 1)
    terminate(1);

  initialize_machine(&argc, &argv);
#ifdef HAVE_QDP
  QDP_initialize(&argc, &argv);
#endif
  g_sync();
  prompt = setup();

  // Load input and run (loop removed)
  if (readin(prompt) != 0) {
    node0_printf("ERROR in readin, aborting\n");
    terminate(1);
  }

  // Set up loop tables
#ifdef IMP
  make_loop_table2();
#endif

  // Start the clocks
  dtime = -dclock();
#ifdef TIMING
  time_dcongrad = 0.0;
  time_fermion_force = 0.0;
  time_fermion_rep = 0.0;
  time_block_nhyp = 0.0;
  time_compute_fhb = 0.0;
//  time_jacobi = 0.0;
#endif

  // Check: compute initial plaquette
  d_plaquette(&ss_plaq, &st_plaq);
  node0_printf("START %.8g %.8g %.8g\n", ss_plaq, st_plaq, ss_plaq + st_plaq);

  f_eps0 = traj_length / (Real)nsteps[0];
  f_eps1 = f_eps0 / (2 * (Real)nsteps[1]);
  g_eps = f_eps1 / (2 * (Real)nsteps[MAX_MASSES]);
  node0_printf("f_eps0 %.4g f_eps1 %.4g g_eps %.4g\n", f_eps0, f_eps1, g_eps);

  // Perform warmup trajectories
  for (traj_done = 0; traj_done < warms; traj_done++)
    update();
  node0_printf("WARMUPS COMPLETED\n");

  // Perform trajectories, reunitarizations and measurements
  for (traj_done = 0; traj_done < trajecs; traj_done++) {
    s_iters = update();
    avs_iters += s_iters;

    /* Do "local" measurements every trajectory! */
    /* The action from the RG trans */
#ifdef IMP  /* For improved action only
               Plaquette action trivially follows from ds[s|t]plaq */
    gauge_action(&ss_plaq);
    node0_printf("ACTION_V %.8g %.8g\n",
                 ss_plaq, ss_plaq / (double)(6.0 * volume));
#endif
    // Measure Polyakov loop and plaquette
    plp = ploop();
    d_plaquette(&ss_plaq, &st_plaq);
    d_plaquette_frep(&ss_plaq_frep, &st_plaq_frep);
#ifdef LU
    // Generate a pseudofermion configuration
    m_iters = f_measure_cl();
#endif
    avm_iters += m_iters;

    node0_printf("GMES %.8g %.8g %d %d %.8g %.8g %.8g %.8g\n",
                 plp.real, plp.imag, s_iters, m_iters,
                 ss_plaq, st_plaq, ss_plaq_frep, st_plaq_frep);

    // Measure every "propinterval" trajectories
    if ((traj_done % propinterval) == (propinterval - 1)) {
      Nmeas++;
#ifdef LU
  #ifdef PCAC
      // Correlators for PCAC relation
      // t direction: Double the lattice via PBC+APBC
      m_iters += pcac_t();
      // x direction: PBC and PBC+APBC
      if (nt < nx)
        m_iters += pcac_x();
  #endif
#endif

      fixflag = NO_GAUGE_FIX;
#ifdef SPECTRUM
  #ifdef SCREEN
      gaugefix(ZUP, 1.5, 500, (Real)GAUGE_FIX_TOL,
               F_OFFSET(staple), F_OFFSET(tempmat1),
               0, NULL, NULL, 0, NULL, NULL);
      fixflag = COULOMB_GAUGE_FIX;
      spect_iters = s_props_cl();
      avspect_iters += spect_iters;
  #else // Spectrum in time direction
      gaugefix(TUP, 1.5, 500, (Real)GAUGE_FIX_TOL,
               F_OFFSET(staple),F_OFFSET(tempmat1),
               0, NULL, NULL, 0, NULL, NULL);
      fixflag = COULOMB_GAUGE_FIX;
      spect_iters = w_spectrum_cl();
      avspect_iters += spect_iters;
  #endif
#endif

    }
    fflush(stdout);
  }
  node0_printf("RUNNING COMPLETED\n");

  // Check: compute final plaquette
  d_plaquette(&ss_plaq, &st_plaq);
  node0_printf("STOP %.8g %.8g %.8g\n",
               ss_plaq, st_plaq, ss_plaq + st_plaq);

  node0_printf("Average cg iters for steps: %.4g\n",
               (double)avs_iters / trajecs);
  if (Nmeas > 0) {
    node0_printf("Average cg iters for measurements: %.4g\n",
                 (double)avm_iters / Nmeas);
#ifdef SPECTRUM
    node0_printf("Average cg iters for spectrum: %.4g\n",
                 (double)avspect_iters / Nmeas);
#endif
  }

  dtime += dclock();
  node0_printf("Time               = %.4g seconds\n", dtime);
#ifdef TIMING
  node0_printf("congrad time       = %.4g seconds\n", time_dcongrad);
  node0_printf("fermion_force time = %.4g seconds\n", time_fermion_force);
//  node0_printf("fermion_rep time   = %.4g seconds\n", time_fermion_rep);
  node0_printf("block_nhyp time    = %.4g seconds\n", time_block_nhyp);
  node0_printf("compute_fhb time   = %.4g seconds\n", time_compute_fhb);
//  node0_printf("jacobi time        = %.4g seconds\n\n", time_jacobi);
#endif
  node0_printf("total_iters = %d\n\n", total_iters);
  fflush(stdout);

  // Save lattice if requested
  if (saveflag != FORGET)
    save_lattice(saveflag, savefile, stringLFN);

  normal_exit(0);
  return 0;
}
// -----------------------------------------------------------------
