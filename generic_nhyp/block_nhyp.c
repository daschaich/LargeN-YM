/********** block_nhyp_arb.c ******************************/
/* Reference:
* Hypercubic smeared links for dynamical fermions.
* By Anna Hasenfratz, Roland Hoffmann, Stefan Schaefer.
* JHEP 0705:029,2007. [hep-lat/0702028]
*/

#include "generic_nhyp_includes.h"

/* parts of this code are specific to NCOL=2,3,4
   valid for any NCOL: calculation of Omega, Staple, and Q
   valid for SU(2,3,4) only: calculation of Q^{-1/2}, including compute_fhb()

T.D. attempt to make this competely general
*/



// -----------------------------------------------------------------
void staple_nhyp(int dir, int dir2, su3_matrix_f *lnk1, su3_matrix_f *lnk2,
                 su3_matrix_f *stp) {

    register int i;
    register site *s;
    msg_tag *tag0, *tag1, *tag2;
    su3_matrix_f tmat;

    // dir is the direction of the original link
    // dir2 is the other direction that defines the staple
    // Get blocked_link[dir2] from direction dir */
    tag0 = start_gather_field(lnk2, sizeof(su3_matrix_f), dir,
                              EVENANDODD, gen_pt[0]);

    // Get blocked_link[dir] from direction dir2 */
    tag1 = start_gather_field(lnk1, sizeof(su3_matrix_f), dir2,
                              EVENANDODD, gen_pt[1]);

    // Start working on the lower staple while we wait for the gathers
    // The lower staple is prepared at x-dir2, stored in tempmat_nhyp1
    // and then gathered to x
    FORALLSITES(i, s)
      mult_su3_an_f(lnk2 + i, lnk1 + i, tempmatf + i);

    // Finish and gather lower staple from direction -dir2
    wait_gather(tag0);
    wait_gather(tag1);
    FORALLSITES(i, s) {
      mult_su3_nn_f(tempmatf + i, (su3_matrix_f *)gen_pt[0][i], &tmat);
      su3mat_copy_f(&tmat, tempmatf + i);
    }

    tag2 = start_gather_field(tempmatf, sizeof(su3_matrix_f),
                              OPP_DIR(dir2), EVENANDODD, gen_pt[2]);

    // Calculate and add upper staple while gather runs
    FORALLSITES(i, s) {
      mult_su3_nn_f(lnk2 + i, (su3_matrix_f *)gen_pt[1][i], &tmat);
      mult_su3_na_sum_f(&tmat, (su3_matrix_f *)gen_pt[0][i], stp + i);
    }

    // Finally add the lower staple
    wait_gather(tag2);
    FORALLSITES(i, s)
      add_su3_matrix_f(stp + i, (su3_matrix_f *)gen_pt[2][i], stp + i);

    cleanup_gather(tag0);
    cleanup_gather(tag1);
    cleanup_gather(tag2);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
#if SMEAR_LEVEL == 3
void block_nhyp1() {
  register int dir, dir2, i;
  register site *st;
  Real f[NCOL], tr1, tr2;
  su3_matrix_f Omega, Q[NCOL];
#if NCOL > 2
  su3_matrix_f eQ;
  int j;
#endif

  tr2 = 1.0 - alpha_smear[2];
  tr1 = alpha_smear[2] / (2.0 * tr2);

  // dir is the direction of the original link
  // dir2 is the other direction that defines the staple
  FORALLUPDIR(dir) {
    FORALLUPDIR(dir2) {
      if (dir2 == dir)
        continue;
      FORALLSITES(i, st)
        clear_su3mat_f(Staple1[dir2][dir] + i);

      // Compute the staple
      staple_nhyp(dir, dir2, gauge_field_thin[dir],
                  gauge_field_thin[dir2], Staple1[dir2][dir]);

      // Make Omega, including IR regulator
      FORALLSITES(i, st) {
        scalar_mult_add_su3_matrix_f(gauge_field_thin[dir] + i,
            Staple1[dir2][dir] + i, tr1, &Q[1]);
        scalar_mult_su3_matrix_f(&Q[1], tr2, &Omega);
        Staple1[dir2][dir][i]=Omega;

        mult_su3_an_f(&Omega, &Omega, &Q[1]);
        scalar_add_diag_su3_f(&Q[1],IR_STAB);
#ifndef NHYP_DEBUG
        compute_fhb(Q[1], f, NULL, 0);
#else
        compute_fhb(&Omega, &Q[1], f, NULL, 0);
#endif

#if NCOL == 2
        scalar_mult_su3_matrix_f(&Omega, f[0], hyplink1[dir2][dir] + i);
#else
        // Compute Q^(-1 / 2)
        scalar_mult_su3_matrix_f(&Q[1], f[1], &eQ);
        for (j = 2; j < NCOL; j++) {
          mult_su3_nn_f(&Q[1], &Q[j-1], &Q[j]);
          scalar_mult_add_su3_matrix_f(&eQ, &Q[j], f[j], &eQ);
        }
        scalar_add_diag_su3_f(&eQ, f[0]);

        // Multiply Omega by eQ = (Omega^dag Omega)^(-1 / 2)
        mult_su3_nn_f(&Omega, &eQ, hyplink1[dir2][dir] + i);
#endif
      }
    }
  }
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
#if (SMEAR_LEVEL>1)
void block_nhyp2() {
  register int dir, dir2, dir3, dir4, i;
  register site *st;
  Real f[NCOL], tr1, tr2;
  su3_matrix_f Omega, Q[NCOL];
#if (NCOL>2)
  su3_matrix_f eQ;
  int j;
#endif

  tr2 = 1.0 - alpha_smear[1];
  tr1 = alpha_smear[1] / (4.0 * tr2);

  FORALLUPDIR(dir) {
    FORALLUPDIR(dir2) {
      if (dir2 == dir)
        continue;

      FORALLSITES(i, st)
        clear_su3mat_f(Staple2[dir2][dir] + i);

      FORALLUPDIR(dir3) {
        if (dir3 == dir || dir3 == dir2)
          continue;
        FORALLUPDIR(dir4) {
          if (dir4 != dir && dir4 != dir2 && dir4 != dir3)
            break;
        }

        // Compute the staple
#if SMEAR_LEVEL == 3
        staple_nhyp(dir, dir3, hyplink1[dir4][dir],
                    hyplink1[dir4][dir3],Staple2[dir2][dir]);
#else // SMEAR_LEVEL == 2
        staple_nhyp(dir, dir3, gauge_field_thin[dir],
                    gauge_field_thin[dir3],Staple2[dir2][dir]);
#endif
      }

      // Make Omega, including IR regulator
      FORALLSITES(i, st) {
        scalar_mult_add_su3_matrix_f(gauge_field_thin[dir] + i,
            Staple2[dir2][dir] + i, tr1, &Q[1]);
        scalar_mult_su3_matrix_f(&Q[1], tr2, &Omega);
        Staple2[dir2][dir][i] = Omega;

        mult_su3_an_f(&Omega, &Omega, &Q[1]);
        scalar_add_diag_su3_f(&Q[1], IR_STAB);
#ifndef NHYP_DEBUG
        compute_fhb(&Q[1], f, NULL, 0);
#else
        compute_fhb(&Omega, &Q[1], f, NULL, 0);
#endif

#if NCOL == 2
        scalar_mult_su3_matrix_f(&Omega, f[0], hyplink2[dir2][dir] + i);
#else
        // Compute Q^(-1 / 2)
        scalar_mult_su3_matrix_f(&Q[1], f[1], &eQ);
        for (j = 2; j < NCOL; j++) {
          mult_su3_nn_f(&Q[1], &Q[j-1], &Q[j]);
          scalar_mult_add_su3_matrix_f(&eQ, &Q[j], f[j], &eQ);
        }
        scalar_add_diag_su3_f(&eQ, f[0]);

        // Multiply Omega by eQ = (Omega^dag Omega)^(-1 / 2)
        mult_su3_nn_f(&Omega, &eQ, hyplink2[dir2][dir] + i);
#endif
      }
    }
  }
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void block_nhyp3() {
  register int dir, dir2, i;
  register site *st;
  Real f[NCOL], tr1, tr2;
  su3_matrix_f Omega, Q[NCOL];
#if NCOL > 2
  su3_matrix_f eQ;
  int j;
#endif

  tr2 = 1.0 - alpha_smear[0];
  tr1 = alpha_smear[0] / (6.0 * tr2);

  FORALLUPDIR(dir) {
    FORALLSITES(i, st)
      clear_su3mat_f(&Staple3[dir][i]);

    FORALLUPDIR(dir2) {
      if (dir2 == dir)
        continue;

      // Compute the staple
#if SMEAR_LEVEL > 1
      staple_nhyp(dir, dir2, hyplink2[dir2][dir], hyplink2[dir][dir2],
                  Staple3[dir]);
#else // SMEAR_LEVEL == 1
      staple_nhyp(dir, dir2, gauge_field_thin[dir], gauge_field_thin[dir2],
                  Staple3[dir]);
#endif
    }

    // Make Omega, including IR regulator
    FORALLSITES(i, st) {
      scalar_mult_add_su3_matrix_f(gauge_field_thin[dir] + i,
                                   Staple3[dir] + i, tr1, &Q[1]);
      scalar_mult_su3_matrix_f(&Q[1], tr2, &Omega);
      Staple3[dir][i] = Omega;
      mult_su3_an_f(&Omega, &Omega, &Q[1]);
      scalar_add_diag_su3_f(&Q[1], IR_STAB);
#ifndef NHYP_DEBUG
      compute_fhb(&Q[1], f, NULL, 0);
#else
      compute_fhb(&Omega, &Q[1], f, NULL, 0);
#endif

#if NCOL == 2
      scalar_mult_su3_matrix_f(&Omega, f[0], gauge_field[dir] + i);
#else
      // Compute Q^(-1 / 2)
      scalar_mult_su3_matrix_f(&Q[1], f[1], &eQ);
      for (j = 2; j < NCOL; j++) {
        mult_su3_nn_f(&Q[1], &Q[j - 1], &Q[j]);
        scalar_mult_add_su3_matrix_f(&eQ, &Q[j], f[j], &eQ);
      }
      scalar_add_diag_su3_f(&eQ, f[0]);

      // Multiply Omega by eQ = (Omega^dag Omega)^(-1 / 2)
      mult_su3_nn_f(&Omega, &eQ, gauge_field[dir] + i);
#endif
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Do SMEAR_LEVEL <= 3 smearing levels, where SMEAR_LEVEL is set in defines.h
void block_nhyp() {
#ifdef TIMING
  TIC(3)
#endif

#if SMEAR_LEVEL == 3
  block_nhyp1();
#endif
#if SMEAR_LEVEL > 1
  block_nhyp2();
#endif
  block_nhyp3();

#ifdef TIMING
  TOC(3, time_block_nhyp)
#endif
}
// -----------------------------------------------------------------
