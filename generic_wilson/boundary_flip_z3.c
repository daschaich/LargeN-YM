/********** boundary_flip_z3.c *************/
/* MIMD version 7 */

/* Multiply timelike links on the last time slice by exp((+/-)2 pi i/NCOL).
   Works on linkf, not link!
   Argument "sign" is PLUS or MINUS, sign of exponent */

/* warning:
does NOT keep track of flips with current_boundary flag */

#include "generic_wilson_includes.h"

void boundary_flip_z3( int sign ) {
  register int i,j,k;
  register site *s;
  register complex phase,tmp;

  phase.real = cos(2*PI/NCOL);
  phase.imag = sin(2*PI/NCOL);
  if(sign == MINUS) phase.imag = - phase.imag;

  FORALLSITES(i,s){
    if(s->t != nt-1)continue;	/* hit only last time slice */
    for(j=0;j<NCOL;j++)for(k=0;k<NCOL;k++){
      CMUL(s->linkf[TUP].e[j][k],phase,tmp);
      s->linkf[TUP].e[j][k] = tmp;
    }
  }
}


