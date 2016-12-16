/********** fermion_rep.c **********************************/
/* MIMD version 7 */

/*
Wrapper for creating link field 'link' in rep of fermions from field 'linkf'
in fundamental rep.  Handles boundary condition stuff,
calls irrep specific routines to perform the actual translation.
*/

#include "generic_wilson_includes.h"

void fermion_rep()  {

   site *st;
   int mu,i;

#ifdef TIMING
TIC(2)
#endif
#ifdef NHYP
/*  block_nhyp looks for the think link in gauge_field_thin[mu],
    and puts the fat link in gauge_field[mu].
    For SF, fixed spatial links at t=0 handled in sf_make_boundary.c
*/
   for(mu=0;mu<4;mu++){
      FORALLDYNLINKS(i,st,mu){
         su3mat_copy_f(&(st->linkf[mu]), gauge_field_thin[mu]+i);
      }
   }

   block_nhyp();
#endif

   for(mu=XUP;mu<=TUP;mu++){
      FORALLDYNLINKS(i,st,mu){
#ifndef NHYP
	 make_fermion_rep_matrix( &(st->linkf[mu]), &(st->link[mu]) );
#else
         make_fermion_rep_matrix( gauge_field[mu]+i, &(st->link[mu]) );
#endif
      }
   }
#ifdef SF
   sf_phases();
#else
/* APBC in time direction     */
   current_boundary=PLUS;
   current_boundary_x=PLUS;
   boundary_flip(MINUS);
#endif

#ifdef TIMING
TOC(2,time_fermion_rep)
#endif

}
