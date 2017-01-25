/* This version uses gathers to get the neighbors */

/* Gauge transform, in place, to axial gauge.
   Then store time-like gauge fields of last time slice in all
   other time slices as well, instead of the unit matrix!
Note that this last step means that this is not really a G.T., it is done
to make the calculation of potentials etc easier */

#include "generic_includes.h"

void ax_gauge() {
  register int i,dir,t,j,k;
  register site *s;
  matrix tmat;
  msg_tag *mtag[4];

  mtag[TUP] = start_gather_site(F_OFFSET(link[TUP]), sizeof(matrix),
      TDOWN, EVENANDODD, gen_pt[TUP] );

  /* Put gauge transformation into staple; it is unity for t=0 */
  FORALLSITES(i,s) if ((s->t)==0 ) {
    for (j=0; j<3; j++)  {
      for (k=0; k<3; k++)  {
        if (j != k)  {
          s->staple.e[j][k] = cmplx(0.0,0.0);
        }
        else  {
          s->staple.e[j][k].real = 1.0;
          s->staple.e[j][k].imag = 0.0;
        }
      }
    }
  }

  /* recursively determine the gauge transformations needed */
  for (t=1;t<nt;t++) {

    wait_gather(mtag[TUP]);

    FORALLSITES(i,s) if ((s->t)==t ) {
      mat_copy((matrix *)(gen_pt[TUP][i]), &(s->staple));
      mat_copy(&(s->link[TUP]), &tmat);
      mult_nn(&(s->staple), &tmat, &(s->link[TUP]));
    }

    if (t<(nt-1)) {
      restart_gather_site(F_OFFSET(link[TUP]), sizeof(matrix),
          TDOWN, EVENANDODD, gen_pt[TUP], mtag[TUP] );
    }
    else{
      cleanup_gather(mtag[TUP]);
      mtag[TUP] = start_gather_site(F_OFFSET(link[TUP]),
          sizeof(matrix), TUP, EVENANDODD, gen_pt[TUP] );
    }
  }

  /* Now do the gauge transformation of the space-like links */
  for (dir=XUP;dir<=ZUP;dir++) {
    mtag[dir] = start_gather_site(F_OFFSET(staple), sizeof(matrix),
        dir, EVENANDODD, gen_pt[dir] );
  }

  for (dir=XUP;dir<=ZUP;dir++) {
    FORALLSITES(i,s) {
      mat_copy(&(s->link[dir]), &tmat);
      mult_nn(&(s->staple), &tmat, &(s->link[dir]));
    }

    wait_gather(mtag[dir]);

    FORALLSITES(i,s) {
      mat_copy(&(s->link[dir]), &tmat);
      mult_na(&tmat, (matrix *)(gen_pt[dir][i]),
          &(s->link[dir]));
    }
  }

  for (dir=XUP;dir<=ZUP;dir++)
    cleanup_gather(mtag[dir]);

  /* recursively copy the last time-like links to all time-slices */
  for (t=nt-2;t>=0;t--) {

    wait_gather(mtag[TUP]);

    FORALLSITES(i,s) if ((s->t)==t ) {
      mat_copy((matrix *)(gen_pt[TUP][i]),&(s->link[TUP]));
    }

    if (t>0) restart_gather_site(F_OFFSET(link[TUP]), sizeof(matrix),
        TUP, EVENANDODD, gen_pt[TUP], mtag[TUP] );
  }

  cleanup_gather(mtag[TUP]);
}
