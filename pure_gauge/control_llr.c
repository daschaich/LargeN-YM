// -----------------------------------------------------------------
// Main procedure for pure-gauge evolution
// Removed hybrid Monte Carlo updates,
// keeping just over-relaxed quasi-heat bath
#define CONTROL
#include "pg_includes.h"


double action();
double actionref(double betaref);
void findEint(double Eint, double delta);
void findEintsmooth(double Eint, double delta);

int main(int argc, char *argv[]) {
  int traj_done, Nint;//, Nmeas = 0;
  int prompt;
  double ss_plaq, st_plaq, dtime, energy;
  complex plp = cmplx(99.0, 99.0);

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

  // Check: compute initial plaquette
  plaquette(&ss_plaq, &st_plaq);
  node0_printf("START %.8g %.8g %.8g\n", ss_plaq, st_plaq, ss_plaq + st_plaq);

  // Perform warmup trajectories
  node0_printf("index 0 %.8g\n", lattice[0].linkf[1].e[0][0].real);
  energy = action();
  node0_printf("GMES %.8g %.8g %.8g %.8g %.8g\n",
               plp.real, plp.imag, ss_plaq, st_plaq, energy);
  //findEintsmooth(Eint,delta);
  node0_printf("index 0 %.8g\n", lattice[0].linkf[1].e[0][0].real);

  srand(iseed);

  // Number of energy intervals, include interval starting at Emax
  Nint = (int)((Emax - Emin) / delta) + 1;
  int Eint, jcount, acounter, k;
  int RMcount;      // Count Robbins--Monro iterations
  double x0[Nint];  // Lower end of energy interval
  double a[Nint];
  double a_i[Njacknife];
  double a_i_new;
  double Reweightexpect; // Reweighted expectation value of the energy
  double *measurement = malloc(trajecs * sizeof(double));
//  double varianz;
  bool Einterval = false;

  for(Eint=0;Eint<Nint;Eint++)
  {
    x0[Eint] = (double)(Eint)*delta + Emin;
    a[Eint] = 0.0;
    for(jcount = 0;jcount<Njacknife;jcount++)
    {
      a_i[jcount] = 1.0;
      a_i_new = 1.0;

      RMcount = 0;
      coldlat();
      Einterval = false;

      for(acounter=0;acounter<ait;acounter++)
      {
        a_i[jcount] = a_i_new;

        if(Einterval == false)
        {
          findEintsmooth(x0[Eint],delta);
          Einterval = true;

          for (traj_done = 0; traj_done < (2*warms); traj_done++)
            updateconst_e(x0[Eint], 1.0);
        }

        for (traj_done = 0; traj_done < warms; traj_done++)
          updateconst_e(x0[Eint], 1.0);

        Reweightexpect=0;
        for(k = 0;k<trajecs;k++)
        {
          measurement[k] = action();
          Reweightexpect += measurement[k];
          updateconst_e(x0[Eint], a_i[jcount]);
        }
        Reweightexpect /= trajecs;
        /*
        // TODO: Should be able to compute this
        //       just by monitoring both <O> and <O^2>
        varianz=0;
        for(k=0;k<trajecs;k++)
        {

        varianz = varianz + pow(measurement[k]-Reweightexpect,2);
        }
        varianz = varianz/trajecs;
        */

        Reweightexpect -= x0[Eint] + 0.5*delta;

        if(RMcount<100)
        {
          a_i_new = a_i[jcount] + 12/(delta*delta)*Reweightexpect;
          node0_printf("a = %.4g off\n", a_i_new);
        }
        else
        {
          a_i_new = a_i[jcount] + 12/(delta*delta*(RMcount+1))*Reweightexpect;
          node0_printf("a = %.4g\n", a_i_new);
        }
        RMcount++;
      }
      a[Eint] += a_i_new / Njacknife;
    }
    //node0_printf("a = %.4g \n", a[Eint]);
    //printf("%d \n",1.0);
  }
  node0_printf("RUNNING COMPLETED\n");

  // Check: compute final plaquette
  plaquette(&ss_plaq, &st_plaq);
  node0_printf("STOP %.8g %.8g %.8g\n", ss_plaq, st_plaq, ss_plaq + st_plaq);

  dtime += dclock();
  node0_printf("Time = %.4g seconds\n", dtime);
  fflush(stdout);

  // Save lattice if requested
  if (saveflag != FORGET)
    save_lattice(saveflag, savefile, stringLFN);

  free(measurement);
  normal_exit(0);
  g_sync();         // Needed by at least some clusters
  return 0;
}


double action() {
  double ssplaq, stplaq, g_act;

  plaquette(&ssplaq, &stplaq);
  g_act = -beta * volume * (ssplaq + stplaq);

  return g_act;
}

double actionref(double betaref) {
  double ssplaq, stplaq, g_act;

  plaquette(&ssplaq, &stplaq);
  g_act = -betaref * volume * (ssplaq + stplaq);

  return g_act;
}


void findEintsmooth(double Eint, double delta) {
  bool Efound = false;
  double energy;
  double energyref;
  double betaref = beta;
  energy = action();
  int counter = 0;
  node0_printf("Eint = %.4g \n", Eint);
  if(energy>=Eint && energy <= (Eint+delta))
  {
    Efound = true;
    beta = betaref;
  }

  while(Efound == false && counter<200)
  {
    update();
    //energy = action();
    energyref = actionref(betaref);
    //node0_printf("energy %.8g %.8g %.8g\n",
      //           energyref, beta, counter);
    //node0_printf("Eint = %.4g \n", Eint);
    //node0_printf("energy %.8g %.8g %.8g\n", energy, Eint, counter);
    //node0_printf("energyref %.8g %.8g %.8g\n", energyref, beta, counter);
    if(energyref>=Eint && energyref <= (Eint+delta))
    {
      Efound = true;
      beta = betaref;
      node0_printf("energyref %.8g %.8g %d\n", energyref, beta, counter);
    }
    else if(energyref>(Eint+delta))
      beta += 0.1;
    else
      beta -= 0.1;

    counter++;
  }
  //beta = betaref;
}





void findEint(double Eint, double delta) {
  bool Efound = false;
  double energy;
  int i, j, a;
  register site *s;
  check_unitarity();
  if(Efound == false)
  {
    FOREVENSITES(i,s)
    {
      for( j = 0;j<4;j++)
      {
        for(a=0;a <NCOL;a++)
          lattice[i].linkf[j].e[a][a] = cmplx(-1.0, 0.0);;

        energy = action();
        node0_printf("energy %.8g\n", energy);
        if(energy>=Eint && energy <= (Eint+delta))
        {
          a=NCOL;
          j=4;
          i=even_sites_on_node;
          Efound = true;
        }
      }
    }
  }
}
// -----------------------------------------------------------------
