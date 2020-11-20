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
  int traj_done;//, Nmeas = 0;
  int prompt;
  double ss_plaq, st_plaq, dtime, energy;
  complex plp = cmplx(99.0, 99.0);
  //double Eint = -3000.0;
  
  double doubleType;
  
  double plp_abs = 0;
  
  //double Emin = -140.0;
  //double Emax = -140.0;
  //double delta=5.0;
  //int Njacknife = 1;
  //int ait=400;
  int seed = iseed;
  
  
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
  node0_printf("index 0 %.8g\n",
                  lattice[0].linkf[1].e[0][0].real);
  energy = action();
    node0_printf("GMES %.8g %.8g %.8g %.8g %.8g\n",
                 plp.real, plp.imag, ss_plaq, st_plaq, energy);
  //findEintsmooth(Eint,delta);
  node0_printf("index 0 %.8g\n",
                  lattice[0].linkf[1].e[0][0].real);
                  
                  
  srand(seed);   
  
  double x0[(int)((Emax-Emin)/delta)];  //lower end of energy interval
	double a[(int)((Emax-Emin)/delta)]; 
	double a_i[Njacknife]; 
	double a_i_new;
	int RobMarchcount; //number of iteration of the rob march algorithm
	double Reweightexpect; //Reweighted expectationvalue of the energy 
	double *meassurement; //smh does not work with meassurement[trajecs]
	double varianz;
  bool Einterval = false; 
  
  //FILE *fpavalue;           
                  
  int Eint = 0;
  int jcount;
  int acounter;
  int k;

  meassurement = malloc(trajecs * sizeof(double));
  
  //fpavalue = fopen("scalar/avalue.txt","w");
  
  for(Eint=0;Eint<=((Emax-Emin)/delta);Eint++)
  {
    x0[Eint] = (double)(Eint)*delta + Emin;
    a[Eint] = 0.0;
    for(jcount = 0;jcount<Njacknife;jcount++)
    {
      a_i[jcount] = 1.0;
		  a_i_new = 1.0;
        
      RobMarchcount = 0;
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
          {
            updateconst_e(x0[Eint],delta,1.0);
          }
        }
        
        for (traj_done = 0; traj_done < (warms); traj_done++)
          {
            updateconst_e(x0[Eint],delta,1.0);
          }
          
        Reweightexpect=0;
				varianz=0;
        
        
        for(k = 0;k<trajecs;k++)
				{
          
					meassurement[k] = action();
					Reweightexpect = Reweightexpect + meassurement[k];
          updateconst_e(x0[Eint],delta,a_i[jcount]);
				}
        Reweightexpect = Reweightexpect/trajecs;
        /*
        for(k=0;k<trajecs;k++)
				{
          
					varianz = varianz + pow(meassurement[k]-Reweightexpect,2);
				}
				varianz = varianz/trajecs;
        */
        
        Reweightexpect = Reweightexpect - x0[Eint] - 0.5*delta;
        
        if(RobMarchcount<100)
				{
					a_i_new = a_i[jcount] + 12/(delta*delta)*Reweightexpect;
           node0_printf("a = %.4g off \n", a_i_new);
				}
				else
				{
					a_i_new = a_i[jcount] + 12/(delta*delta*(RobMarchcount+1))*Reweightexpect;
           node0_printf("a = %.4g \n", a_i_new);
				}
        RobMarchcount = RobMarchcount + 1;
        
      }
      
      a[Eint] = a[Eint]+a_i_new/Njacknife;
      
    }
    //node0_printf("a = %.4g \n", a[Eint]);
    //printf("%d \n",1.0);
    //fprintf(fpavalue,"Enter a sentence:\n");
    //fprintf(fpavalue, "%d \n",1.0);
  }
  //fclose(fpavalue);
  /*
  for (traj_done = 0; traj_done < warms; traj_done++)
  {
    update();
    energy = action();
    node0_printf("GMES %.8g %.8g %.8g %.8g %.8g\n",
                 plp.real, plp.imag, ss_plaq, st_plaq, energy);
  }
  node0_printf("WARMUPS COMPLETED\n");

  // Perform trajectories, reunitarizations and measurements
  for (traj_done = 0; traj_done < trajecs; traj_done++) {
    //updateconst_e(-3000.0,delta);
    //update();
    energy = action();
    
    if (traj_done > 499)
    {
      plp_abs = plp_abs + sqrt(plp.real*plp.real + plp.imag*plp.imag);
    }
    
    // Measure and print Polyakov loop and plaquette
    // after every trajectory
    plaquette(&ss_plaq, &st_plaq);
    plp = ploop(TUP);
    node0_printf("GMES %.8g %.8g %.8g %.8g %.8g\n",
                 plp.real, plp.imag, ss_plaq, st_plaq, energy);
                 
    
    fflush(stdout);

    // More expensive measurements every "propinterval" trajectories
    // Measure every "propinterval" trajectories
    if ((traj_done % propinterval) == (propinterval - 1)) {
//      Nmeas++;
      // Nothing yet...
    }
  }
  plp_abs = plp_abs/500;
  node0_printf("Polabsolute %.8g\n", plp_abs);
  */
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

  free(meassurement);
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
      node0_printf("energyref %.8g %.8g %.8g\n", energyref, beta, counter);
    }
    else if(energyref>(Eint+delta))
    {
      beta = beta+0.1;
    }
    else
    {
      beta = beta-0.1;
    }
    counter = counter +1;
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
    //for( i = 0;i<volume;i=i+2)
    //FOREVENSITES(i,s)
    //for( i = 0;i<even_sites_on_node;i=i+1)
    for( i = 0;i<volume;i=i+2)
    {
      for( j = 0;j<4;j++)
      {
        for(a=0;a <NCOL;a++)
        {
          
          lattice[i].linkf[j].e[a][a] = cmplx(-1.0, 0.0);;
          
        }
        energy = action();
          node0_printf("energy %.8g\n",
                 energy);
          if(energy>=Eint && energy <= (Eint+delta))
          {
            a=NCOL;
            j=4;
            i=volume;
            //i=even_sites_on_node;
            Efound = true;
          }
      }
    }
  }
}
// -----------------------------------------------------------------
