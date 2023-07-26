// -----------------------------------------------------------------
// Main procedure and helpers for pure-gauge over-relaxed quasi-heatbath
#define CONTROL
#include "pg_includes.h"

int main(int argc, char *argv[]) {
  int prompt;
  int traj_done, RMcount, Ncount, Intcount;//, Nmeas = 0;
  double ss_plaq, st_plaq, E, dtime, save_a, rate;
  double Reweightexpect;    // Reweighted expectation value of the energy
  double ReweightSquareexpect;
  double variance;
  double CurrentE;
  

  // Set up
  setlinebuf(stdout); // DEBUG
  initialize_machine(&argc,&argv);
  g_sync();
  if (remap_stdio_from_args(argc, argv) == 1)
    terminate(1);

  // Load input and run (loop removed)
  prompt = setup();
  if (readin(prompt) != 0) {
    node0_printf("ERROR in readin, aborting\n");
    terminate(1);
  }
  dtime = -dclock();
  
  int nrintervals = (int)((Emax-Emin)/delta)+1;
  node0_printf("nrintervals %.8g \n",
               Emax);
  double Eint[nrintervals];
  double aint[nrintervals];
  double accrate_it = 0;
  // Monitor overall acceptance in monteconst_e.c
  accept = 0;
  reject = 0;

  // Check: compute initial plaquette and energy
  plaquette(&ss_plaq, &st_plaq);
  E = gauge_action();
  node0_printf("START %.8g %.8g %.8g %.8g\n",
               ss_plaq, st_plaq, ss_plaq + st_plaq, E);
  save_a = a;
  for(Intcount = 0; Intcount <= (int)((Emax-Emin)/delta); Intcount++) {
    aint[Intcount] = 0;
    Eint[Intcount] = Emin + Intcount*delta;
    for(Ncount = 0; Ncount < Njacknife; Ncount++) {
      a = 1.0;
      constrained = 0;
      startlat_p = reload_lattice(startflag, startfile);
      // Check: compute initial plaquette and energy
      plaquette(&ss_plaq, &st_plaq);
      E = gauge_action();
      node0_printf("START %.8g %.8g %.8g %.8g\n",
               ss_plaq, st_plaq, ss_plaq + st_plaq, E);

      // Unconstrained warmup sweeps before searching for energy interval
      //for (traj_done = 0; traj_done < 10; traj_done++)
      //  update();
      //node0_printf("WARMUPS COMPLETED\n");

      // Terminates if interval not found
      constrained = 0;
      a = 1.0;
      findEint(Eint[Intcount]+0.5*delta-0.5/1.0*delta);
      //a = save_a;
      // Robbins--Monro (RM) iterations
      constrained = 1;
      for (RMcount = 0; RMcount < ait; RMcount++) {
        // Constrained warm-up sweeps in each RM iteration, with a=1
        //save_a = a;
        //a = 1.0;
        accrate_it=0.0;
        accept = 0;
        reject = 0;
        for (traj_done = 0; traj_done < warms; traj_done++)
          updateconst_e(Eint[Intcount]);
        //a = save_a;

        Reweightexpect = 0.0;
        ReweightSquareexpect = 0.0;
        for (traj_done = 0; traj_done < trajecs; traj_done++) {

          updateconst_e(Eint[Intcount]);
          // Accumulate after update
          CurrentE = gauge_action()-Eint[Intcount] - 0.5 * delta;
          Reweightexpect += CurrentE;
          ReweightSquareexpect += CurrentE*CurrentE;

          // More expensive measurements every "measinterval" sweeps
          //if ((traj_done % measinterval) == (measinterval - 1)) {
//            Nmeas++;
            //   Nothing yet...
            //}
        }
        Reweightexpect /= trajecs;
        ReweightSquareexpect /= trajecs;
        //ReweightSquareexpect = ReweightSquareexpect - 2.0*(Eint[Intcount]+delta/2.0)*Reweightexpect + Eint[Intcount]*Eint[Intcount] + Eint[Intcount]*delta + delta*delta/4.0; 
        //Reweightexpect -= Eint[Intcount] + 0.5 * delta;
        
        variance = ReweightSquareexpect - Reweightexpect*Reweightexpect;
        
        accrate_it = ((double)accept/((double)(accept + reject)));
        node0_printf("Acc rate this step %.8g\n", accrate_it);
        
        if(((double)accept/((double)(accept + reject)))<0.5){
          hmc_steps = hmc_steps + 20;
          RMcount = RMcount - 1;
          node0_printf("Nr hmc steps changed to %d\n", hmc_steps);
          node0_printf("redo step");
        }
        else if(((double)accept/((double)(accept + reject)))>0.9){
          if( (hmc_steps - 10)>0 ){
            hmc_steps = hmc_steps - 10;
            node0_printf("Nr hmc steps changed to %d\n", hmc_steps);
          }
        }
        
          
        
        // Hard-code under-relaxation to begin after 25 RM iterations
        if(accrate_it>=0.5){
          if (RMcount < 10){
              if(variance>deltaSq){
                a += 1.0*1.0 * Reweightexpect / (variance);
              }
              else{
                a += 1.0*12.0 * Reweightexpect / (deltaSq);
              }
            
          }
          else{
            
              if(variance>deltaSq){
                a += 1.0*1.0 * Reweightexpect / (variance * (RMcount - 10.0 + 1.0));
              }
              else{
                a += 1.0*12.0 * Reweightexpect / (deltaSq * (RMcount - 10.0 + 1.0));
              }
            
          }
        }
        node0_printf("RM ITER %d Reweightexpect %.8g a %.8g\n", RMcount + 1, Reweightexpect,a);
        // TODO: I think acceptance rate for each RM iteration
        //       would be more interesting than the overall one below...
      }
      aint[Intcount] = aint[Intcount] + a/((double)(Njacknife));
    }
  }
  node0_printf("RUNNING COMPLETED\n");

  // Check: compute final plaquette and energy
  plaquette(&ss_plaq, &st_plaq);
  E = gauge_action();
  node0_printf("STOP %.8g %.8g %.8g %.8g\n",
               ss_plaq, st_plaq, ss_plaq + st_plaq, E);

  //rate = (double)accept / ((double)(accept + reject));
  //node0_printf("Overall acceptance %d of %d = %.4g\n",
               //accept, accept + reject, rate);
  dtime += dclock();
  node0_printf("Time = %.4g seconds\n", dtime);
  fflush(stdout);

  // Save lattice if requested
  if (saveflag != FORGET)
    save_lattice(saveflag, savefile, stringLFN);

  normal_exit(0);         // Needed by at least some clusters
  return 0;
}
// -----------------------------------------------------------------
