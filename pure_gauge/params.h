// -----------------------------------------------------------------
// Structure for passing simulation parameters to each node
#ifndef _PARAMS_H
#define _PARAMS_H
#include "../include/macros.h"  // For MAXFILENAME

typedef struct {
  int stopflag;           // 1 if it is time to stop
  int nx, ny, nz, nt;     // Lattice dimensions
  int iseed;  /* for random numbers */

  int warms;  /* the number of warmup trajectories */
  int trajecs;  /* the number of real trajectories */
  int steps;  /* number of steps for updating */
  int stepsQ; /* number of steps for qhb */
  int propinterval;     /* number of trajectories between measurements */
  int startflag;  /* what to do for beginning lattice */
  int saveflag;   /* what to do with lattice at end */
  Real beta;  /* gauge coupling */
  Real epsilon; /* time step */
  char startfile[MAXFILENAME], savefile[MAXFILENAME];
  char stringLFN[MAXFILENAME];  /** ILDG LFN if applicable ***/
  double Emin;
  double Emax;
  double delta;
  int ait;
  int Njacknife;
}  params;
#endif
// -----------------------------------------------------------------
