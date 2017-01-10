/* This version uses gathers to get the neighbors */

/* Measure the average plaquette of the space-space and
   space-time plaquettes */

#include "generic_includes.h"

void d_plaquette(double *ss_plaq,double *st_plaq) {
  register int i,dir1,dir2;
  register site *s;
  register su3_matrix *m1,*m4;
  double ss_sum,st_sum;
  msg_tag *mtag0,*mtag1;
  ss_sum = st_sum = 0.0;
  for(dir1=YUP;dir1<=TUP;dir1++){
    for(dir2=XUP;dir2<dir1;dir2++){

      mtag0 = start_gather_site(F_OFFSET(link[dir2]), sizeof(su3_matrix),
          dir1, EVENANDODD, gen_pt[0]);
      mtag1 = start_gather_site(F_OFFSET(link[dir1]), sizeof(su3_matrix),
          dir2, EVENANDODD, gen_pt[1]);

      FORALLSITES(i,s){
        m1 = &(s->link[dir1]);
        m4 = &(s->link[dir2]);
        mult_su3_an(m4,m1,&(s->tempmat1));
      }

      wait_gather(mtag0);
      FORALLSITES(i,s){
        mult_su3_nn(&(s->tempmat1),(su3_matrix *)(gen_pt[0][i]),
            &(s->staple));
      }

      wait_gather(mtag1);
      FORALLSITES(i,s){
        if(dir1==TUP)st_sum += (double)
          realtrace_su3((su3_matrix *)(gen_pt[1][i]),&(s->staple));
        else              ss_sum += (double)
          realtrace_su3((su3_matrix *)(gen_pt[1][i]),&(s->staple));
      }

      cleanup_gather(mtag0);
      cleanup_gather(mtag1);
    }
  }
  g_doublesum(&ss_sum);
  g_doublesum(&st_sum);
  /**{Real xxx; xxx = ss_sum; g_floatsum(&xxx); ss_sum=xxx;}
    {Real xxx; xxx = st_sum; g_floatsum(&xxx); st_sum=xxx;}**/
  *ss_plaq = ss_sum / ((double)(3.0 * volume));
  *st_plaq = st_sum / ((double)(3.0 * volume));
}
