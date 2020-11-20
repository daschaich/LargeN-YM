// -----------------------------------------------------------------
// Structure for passing simulation parameters to each node
#ifndef _PARAMS_H
#define _PARAMS_H
#include "../include/macros.h"  // For MAXFILENAME

typedef struct {
  int stopflag;           // 1 if it is time to stop

  // Initialization parameters
  int nx, ny, nz, nt;     // Lattice dimensions
  int iseed;              // For random numbers

  // Over-relaxation parameters
  int warms;              // The number of warmup trajectories
  int trajecs;            // The number of real trajectories
  int steps;              // Number of over-relaxation steps per trajectory
  int stepsQ;             // Number of quasi-heatbath steps per trajectory
  int propinterval;       // Number of trajectories between measurements
  int startflag;          // What to do for beginning lattice
  int saveflag;           // What to do with lattice at end
  Real beta;              // Gauge coupling
  char startfile[MAXFILENAME], savefile[MAXFILENAME];
  char stringLFN[MAXFILENAME];  // ILDG LFN if applicable

#ifdef LLR
  // LLR parameters
  Real Emin, Emax;    // Energy range to scan
  Real delta;         // Size of energy interval
  int ait;            // Number of Robbins--Monro iterations
  int Njacknife;      // Number of repetitions to jackknife or bootstrap
#endif
}  params;
#endif
// -----------------------------------------------------------------
