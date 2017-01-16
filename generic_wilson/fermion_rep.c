// -----------------------------------------------------------------
// Wrapper for creating DIMFxDIMF link matrices from NCOLxNCOL linkf matrices
// Handles boundary condition stuff
// Calls irrep-specific routines to perform the actual translation.
#include "generic_wilson_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void fermion_rep() {
  register int mu, i;
  register site *s;

#ifdef TIMING
  TIC(2);
#endif
  // block_nhyp looks for the thin link in gauge_field_thin[mu],
  // and puts the smeared link in gauge_field[mu]
  FORALLUPDIR(mu) {
    FORALLSITES(i, s)
      su3mat_copy_f(&(s->linkf[mu]), gauge_field_thin[mu] + i);
  }
  block_nhyp();

  FORALLUPDIR(mu) {
    FORALLSITES(i, s)
      make_fermion_rep_matrix(gauge_field[mu] + i, &(s->link[mu]));
  }
  // Anti-periodic BCs in time direction
  current_boundary = PLUS;
  boundary_flip(MINUS);

#ifdef TIMING
  TOC(2, time_fermion_rep);
#endif
}
// -----------------------------------------------------------------
