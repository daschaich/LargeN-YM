/*** sf_phases *********************************************************
   twisted b.c. in spatial directions for SF;
   imposed by multiplying every link[dir] by the same phase
*/

#include "generic_wilson_includes.h"

void sf_phases(){
site *st;
int mu,is;

   for(mu=XUP;mu<=ZUP;mu++){
      FORALLDYNLINKS(is,st,mu){
            sf_phase( &(st->link[mu]), mu);
      }
   }
}

void sf_phase(su3_matrix *a, int dir){
register int i,j;
register complex temp;

   if(dir>ZUP){
      if(this_node==0)printf("Illegal value %d for spatial dir\n", dir);
      terminate(1);
   }

   for(i=0;i<DIMF;i++)for(j=0;j<DIMF;j++){
      CMUL( a->e[i][j], ferm_phases[dir], temp );
      a->e[i][j]=temp;
   }
}
