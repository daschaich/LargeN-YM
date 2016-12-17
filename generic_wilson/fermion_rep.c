// -----------------------------------------------------------------
// Wrapper for creating DIMFxDIMF matrices 'link' from NCOLxNCOL matrices 'linkf'
// Handles boundary condition stuff
// Calls irrep-specific routines to perform the actual translation.
#include "generic_wilson_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void fermion_rep() {
  register int mu, i;
  register site *st;

#ifdef TIMING
  TIC(2)
#endif
#ifdef NHYP
    /*  block_nhyp looks for the thin link in gauge_field_thin[mu],
        and puts the smeared link in gauge_field[mu].
        */
    for(mu=0;mu<4;mu++){
      FORALLSITES(i,st)
        su3mat_copy_f(&(st->linkf[mu]), gauge_field_thin[mu]+i);
    }

  block_nhyp();
#endif

  for(mu=XUP;mu<=TUP;mu++){
    FORALLSITES(i,st) {
#ifndef NHYP
      make_fermion_rep_matrix(&(st->linkf[mu]), &(st->link[mu]));
#else
      make_fermion_rep_matrix(gauge_field[mu]+i, &(st->link[mu]));
#endif
    }
  }
  // Anti-periodic BCs in time direction
  current_boundary=PLUS;
  current_boundary_x=PLUS;
  boundary_flip(MINUS);

#ifdef TIMING
  TOC(2,time_fermion_rep)
#endif
}
// -----------------------------------------------------------------
