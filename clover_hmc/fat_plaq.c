/*** do measurements with the fat fundamental link ******************/
void fat_plaq() {

site *st;
int mu,i;
double dssplaq,dstplaq;
complex plp;

   /* copy the fat link into st->linkf */
   for(mu=0;mu<4;mu++){
      FORALLDYNLINKS(i,st,mu){
         su3mat_copy_f( gauge_field[mu]+i, &(st->linkf[mu]) );
      }
   }

   /* call the Polyakov loop measuring program */
   plp = ploop();

   /* call plaquette measuring process */
   d_plaquette(&dssplaq,&dstplaq);

   if(this_node==0)printf("GFAT %e %e %e %e\n",
        (double)plp.real,(double)plp.imag,dssplaq,dstplaq);

   /* restore the thin link */
   for(mu=0;mu<4;mu++){
      FORALLDYNLINKS(i,st,mu){
         su3mat_copy_f( gauge_field_thin[mu]+i, &(st->linkf[mu]) );
      }
   }
} 
