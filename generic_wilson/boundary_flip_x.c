/********** boundary_flip_x.c *************/
/* MIMD version 7 */

/* Flip the x direction boundary conditions.
   (Just multiply x links on last x slice by -1.)
   argument "sign" is PLUS or MINUS, says what you want the boundary
   conditions to become.  If the boundary conditions are already set
   to "sign", you get a warning.
*/

#include "generic_wilson_includes.h"

void boundary_flip_x( int sign ) {
  register int i,j,k;
  register site *s;

  if(this_node==0 && sign==current_boundary_x)
    printf("WARNING: you lost track of the x boundary conditions!\n");
  if(sign==current_boundary_x)return;


  FORALLSITES(i,s){
    if(s->x != nx-1)continue;	/* hit only last x slice */
    for(j=0;j<DIMF;j++)for(k=0;k<DIMF;k++){
      s->link[XUP].e[j][k].real *= -1.0;
      s->link[XUP].e[j][k].imag *= -1.0;
    }
  }

  current_boundary_x = sign;
}


