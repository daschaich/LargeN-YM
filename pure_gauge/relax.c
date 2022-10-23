// -----------------------------------------------------------------
// Microcanonical over-relaxation by doing successive SU(2) gauge hits
#include "pg_includes.h"

void relax() {
  register int dir, i;
  register site *s;
  int istep, Nhit, subgrp, ina, inb, count;
  int parity, index_a[N_OFFDIAG], index_b[N_OFFDIAG];
  Real a0, a1, a2, a3, asq, norm;
  su2_matrix u;
  matrix action;

  Nhit = (int)N_OFFDIAG;    // NCOL * (NCOL - 1) / 2
  // Set up SU(2) subgroup indices [a][b], always with a < b
  count = 0;
  for (ina = 0; ina < NCOL; ina++) {
    for (inb = ina + 1; inb < NCOL; inb++) {
      index_a[count] = ina;
      index_b[count] = inb;
      count++;
    }
  }
  if (count != Nhit) {      // Sanity check
    node0_printf("ERROR: %d rather than %d subgroups found", count, Nhit);
    terminate(1);
  }

  for (istep = 0 ; istep < steps; istep++) {
    for (parity = ODD; parity <= EVEN; parity++) {
      FORALLUPDIR(dir) {
        // Compute the gauge force (updating every s->staple)
        dsdu_qhb(dir, parity);

        // Now for the overrelaxed updating
        for (subgrp = 0; subgrp < Nhit; subgrp++) {
          // Pick out this SU(2) subgroup
          ina = index_a[subgrp];
          inb = index_b[subgrp];
          FORSOMEPARITY(i, s, parity) {
            mult_na(&(s->linkf[dir]), &(s->staple), &action);

            // Decompose the action into SU(2) subgroups
            // using Pauli matrix expansion
            // The SU(2) hit matrix is represented as
            //    a0 + i * Sum j (sigma j * aj)
            a0 = action.e[ina][ina].real + action.e[inb][inb].real;
            a3 = action.e[ina][ina].imag - action.e[inb][inb].imag;
            a1 = action.e[ina][inb].imag + action.e[inb][ina].imag;
            a2 = action.e[ina][inb].real - action.e[inb][ina].real;

            // Normalize and complex conjugate u
            asq = a0 * a0 + a1 * a1 + a2 * a2 + a3 * a3;
            norm = 1.0 / sqrt((double)asq);
            a0 *=  norm;
            a1 *= -norm;
            a2 *= -norm;
            a3 *= -norm;

            // Elements of SU(2) matrix
            u.e[0][0] = cmplx( a0, a3);
            u.e[0][1] = cmplx( a2, a1);
            u.e[1][0] = cmplx(-a2, a1);
            u.e[1][1] = cmplx( a0,-a3);

            // Do SU(2) hit on all links twice (to overrelax)
            left_su2_hit_n(&u, ina, inb, &(s->linkf[dir]));
            left_su2_hit_n(&u, ina, inb, &(s->linkf[dir]));
          }
        }
      }
    }
  }
}
// -----------------------------------------------------------------
