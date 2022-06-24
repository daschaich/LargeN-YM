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
  int warms;              // The number of warmup sweeps
  int sweeps;             // The number of sweeps with measurements
  int steps;              // Number of over-relaxation steps per sweep
  int stepsQ;             // Number of quasi-heatbath steps per sweep
  int nsteps;             // Number of hmc steps per trajectory
  Real traj_length;       // Length of trajectory of hmc
  Real eps;
  int startflag;          // What to do for beginning lattice
  int saveflag;           // What to do with lattice at end
  Real beta;              // Gauge coupling
  char startfile[MAXFILENAME], savefile[MAXFILENAME];
  char stringLFN[MAXFILENAME];  // ILDG LFN if applicable

#ifndef LLR
  int measinterval;       // Number of sweepsbetween measurements
#endif

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
