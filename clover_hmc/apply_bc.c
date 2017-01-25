/*
This routine takes care of boundary conditions and twist angles,
as they affect the chain rule.
Applies irrespective of the fermion representation
*/

#include "cl_dyn_includes.h"

void apply_bc(matrix_f *sigmaf, int dir, int t){

  int i,j;
  if( t==nt-1 && dir==TUP && current_boundary==MINUS ){
    for(i=0;i<NCOL;i++)for(j=0;j<NCOL;j++){
      sigmaf->e[i][j].real *= -1.0;
      sigmaf->e[i][j].imag *= -1.0;
    }
  }
}
