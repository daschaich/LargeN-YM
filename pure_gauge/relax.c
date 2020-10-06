// -----------------------------------------------------------------
// Microcanonical over-relaxation by doing successive SU(2) gauge hits
#include "pg_includes.h"

void relax(int NumStp) {
  register int dir, i;
  register site *st;
  int NumTrj, Nhit, subgrp, ina, inb, j, count;
  int parity, index_a[N_OFFDIAG], index_b[N_OFFDIAG];
  Real a0,a1,a2,a3,asq,r;
  su2_matrix u;
  matrix_f action;

  Nhit = (int)N_OFFDIAG;    // NCOL * (NCOL - 1) / 2
  // Set up SU(2) subgroup indices [a][b], always with a < b
  count = 0;
  for (i = 0; i < NCOL; i++) {
    for (j = i + 1; j < NCOL; j++) {
      index_a[count] = i;
      index_b[count] = j;
      count++;
    }
  }
  if (count != Nhit) {      // Sanity check
    node0_printf("ERROR: %d rather than %d subgroups found", count, Nhit);
    terminate(1);
  }

  for( NumTrj = 0 ; NumTrj < NumStp; NumTrj++) {
    for(parity=ODD;parity<=EVEN;parity++) {
      FORALLUPDIR(dir) {
        // Compute the gauge force
        dsdu_qhb(dir, parity);

        // Now for the overrelaxed updating
        for (subgrp = 0; subgrp < Nhit; subgrp++) {
          // Pick out this SU(2) subgroup
          ina = index_a[subgrp];
          inb = index_b[subgrp];
          FORSOMEPARITY(i, st, parity) {
            mult_na_f(&(st->linkf[dir]), &(st->staple), &action);

            /*decompose the action into SU(2) subgroups using Pauli matrix expansion */
            /* The SU(2) hit matrix is represented as a0 + i * Sum j (sigma j * aj)*/
            a0 =  action.e[ina][ina].real + action.e[inb][inb].real;
            a3 =  action.e[ina][ina].imag - action.e[inb][inb].imag;
            a1 =  action.e[ina][inb].imag + action.e[inb][ina].imag;
            a2 =  action.e[ina][inb].real - action.e[inb][ina].real;

            /* Normalize and complex conjugate u */
            asq = a0*a0 + a1*a1 + a2*a2 + a3*a3;
            r = sqrt((double)asq );
            a0 = a0/r; a1 = -a1/r; a2 = -a2/r; a3 = -a3/r;

            /* Elements of SU(2) matrix */
            u.e[0][0] = cmplx( a0, a3);
            u.e[0][1] = cmplx( a2, a1);
            u.e[1][0] = cmplx(-a2, a1);
            u.e[1][1] = cmplx( a0,-a3);

            /* Do SU(2) hit on all links twice (to overrelax)  */
            left_su2_hit_n_f(&u,ina,inb,&(st->linkf[dir]));
            left_su2_hit_n_f(&u,ina,inb,&(st->linkf[dir]));
          } /*   st */
        } /*  hits */
      } /*  direction */
    }} /* parity, NumTrj*/
} /* relax */
// -----------------------------------------------------------------
