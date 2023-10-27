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
  int startflag;          // What to do for beginning lattice
  int saveflag;           // What to do with lattice at end
  char startfile[MAXFILENAME], savefile[MAXFILENAME];
  char stringLFN[MAXFILENAME];  // ILDG LFN if applicable

  int ora_steps;          // Number of over-relaxation steps per sweep
  int qhb_steps;          // Number of quasi-heatbath steps per sweep
#ifdef HMC
  int hmc_steps;          // Number of hmc steps per trajectory
  Real traj_length;       // Trajectory length
#endif

#ifdef LLR
  Real Emin;              // Lower edge of lowest energy interval
  Real Emax;              // Lower edge of highest energy interval
  Real delta;             // Size of energy interval
  Real C_Gauss;           // Coefficient of Gaussian window function
  int NRiter;             // Number of Newton--Raphson iterations
  int RMiter;             // Number of subsequent Robbins--Monro iterations
  int Nj;                 // Number of Jackknife samples
#endif
} params;
#endif
// -----------------------------------------------------------------
