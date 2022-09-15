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

  // Generic update parameters
  int warms;              // The number of warmup trajectories/sweeps
  int trajecs;            // The number of trajectories with measurements
  int measinterval;       // Number of trajectories between measurements
  Real beta;              // Gauge coupling
  Real a;                 // LLR observable, otherwise fixed to a=1
  int startflag;          // What to do for beginning lattice
  int saveflag;           // What to do with lattice at end
  char startfile[MAXFILENAME], savefile[MAXFILENAME];
  char stringLFN[MAXFILENAME];  // ILDG LFN if applicable

#ifndef HMC
  // Over-relaxation parameters
  int ora_steps;              // Number of over-relaxation steps per sweep
  int qhb_steps;              // Number of quasi-heatbath steps per sweep
#else
  // HMC parameters
  int hmc_steps;              // Number of hmc steps per trajectory
  Real traj_length;           // Trajectory length
  Real eps;
#endif

#ifdef LLR
  // LLR parameters
  Real Emin;          // Lower edge of energy interval
  Real delta;         // Size of energy interval
  int ait;            // Number of Robbins--Monro iterations
#endif
}  params;
#endif
// -----------------------------------------------------------------
