// -----------------------------------------------------------------
// Derivative of the nHYP link w.r.t. the constituent links
// Two contributions:
// 1) The (thin) link itself : added to Sigma in Sigma_update1
// 2) The (fat) staples      : made in compute-sigma23 and put in SimgaH
#include "cl_dyn_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Helper routine for Sigma_update1
void make_2hermitian_f(su3_matrix_f *A) {
  int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = i; j < NCOL; j++) {
      A->e[i][j].real =  A->e[i][j].real + A->e[j][i].real;
      A->e[i][j].imag =  A->e[i][j].imag - A->e[j][i].imag;
      A->e[j][i].real =  A->e[i][j].real;
      A->e[j][i].imag = -A->e[i][j].imag;
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute Gamma
void Sigma_update1(const int dir, su3_matrix_f* sigma_off, su3_matrix_f *stp,
                   su3_matrix_f **lambda, const Real alpha1,
                   const Real alpha2, const int sigfresh) {

  register int i;
  register site *s;
  Real f[NCOL], bb[NCOL][NCOL];
  complex traces[NCOL], tc;
  su3_matrix_f Gamma, Qisqrt, Q, Omega, tmat, SigmaOmega;
#if NCOL > 2
  int j;
  complex tc2, tc3;
  su3_matrix_f Q2;
#if NCOL > 3
  complex tc4;
  su3_matrix_f Q3, tmat2;
#endif
#endif

  FORALLSITES(i, s) {
    // Make Omega, Q, Q^2
    Omega = stp[i];
    mult_su3_an_f(&Omega, &Omega, &Q);
    // IR regulator (set in defines.h) doesn't change derivatives of Q
    scalar_add_diag_su3_f(&Q, IR_STAB);

    // Compute inverse sqrt
#ifndef NHYP_DEBUG
    compute_fhb(&Q, f, bb, 1);
#else
    compute_fhb(&Omega, &Q, f, bb, 1);
#endif

#if NCOL == 2
    clear_su3mat_f(&Qisqrt);
    Qisqrt.e[0][0].real = f[0];
    Qisqrt.e[1][1].real = f[0];
#else   // NCOL > 2
    mult_su3_nn_f(&Q, &Q, &Q2);
    scalar_mult_su3_matrix_f(&Q, f[1], &Qisqrt);
    scalar_mult_sum_su3_matrix_f(&Q2, f[2], &Qisqrt);
#if NCOL > 3
    mult_su3_nn_f(&Q, &Q2, &Q3);
    scalar_mult_sum_su3_matrix_f(&Q3, f[3], &Qisqrt);
#endif
    scalar_add_diag_su3_f(&Qisqrt, f[0]);
#endif

    // We'll need Sigma*Omega a few times
    mult_su3_nn_f(sigma_off + i, &Omega, &SigmaOmega);

    /* now the B matrices and their traces with Sigma*Omega*/
    tc=trace_su3_f(&SigmaOmega);
#if (NCOL==2)
    traces[0].real = tc.real * bb[0][0];
    traces[0].imag = tc.imag * bb[0][0];
    clear_su3mat_f(&Gamma);
    c_scalar_add_diag_su3_f(&Gamma, &traces[0]);
#else
    tc2 = complextrace_su3_f(&Q, &SigmaOmega);
    tc3 = complextrace_su3_f(&Q2, &SigmaOmega);
#if NCOL > 3
    tc4 = complextrace_su3_f(&Q3, &SigmaOmega);
#endif
    for (j = 0; j < NCOL; j++) {
      traces[j].real = tc.real * bb[j][0] + tc2.real * bb[j][1]
                                          + tc3.real * bb[j][2];
      traces[j].imag = tc.imag * bb[j][0] + tc2.imag * bb[j][1]
                                          + tc3.imag * bb[j][2];
#if NCOL > 3
      traces[j].real += tc4.real * bb[j][3];
      traces[j].imag += tc4.imag * bb[j][3];
#endif
    }

    // The contributions to A tr(B_i Sigma Omega) Q^(i)
    c_scalar_mult_su3mat_f(&Q, &traces[1], &Gamma);
    c_scalar_mult_sum_su3mat_f(&Q2, &traces[2], &Gamma);
#if NCOL > 3
    c_scalar_mult_sum_su3mat_f(&Q3, &traces[3], &Gamma);
#endif
    c_scalar_add_diag_su3_f(&Gamma, &traces[0]);

    // The terms proportional to f_i
    scalar_mult_sum_su3_matrix_f(&SigmaOmega, f[1], &Gamma);
    mult_su3_nn_f(&SigmaOmega, &Q, &tmat);
    scalar_mult_sum_su3_matrix_f(&tmat, f[2], &Gamma);
    mult_su3_nn_f(&Q, &SigmaOmega, &tmat);
    scalar_mult_sum_su3_matrix_f(&tmat, f[2], &Gamma);
#if NCOL > 3
    mult_su3_nn_f(&tmat, &Q, &tmat2);
    scalar_mult_sum_su3_matrix_f(&tmat2, f[3], &Gamma);
    mult_su3_nn_f(&SigmaOmega, &Q2, &tmat);
    scalar_mult_sum_su3_matrix_f(&tmat, f[3], &Gamma);
    mult_su3_nn_f(&Q2, &SigmaOmega, &tmat);
    scalar_mult_sum_su3_matrix_f(&tmat, f[3], &Gamma);
#endif

#endif
    /* Gamma = (A + Adag)Omega^dag + Q^{-1 / 2}.Sigma */
    make_2hermitian_f(&Gamma);
    mult_su3_na_f(&Gamma, &Omega, &tmat);

    mult_su3_nn_f(&Qisqrt, sigma_off + i, &Gamma);
    sum_su3_matrix_f(&tmat, &Gamma);
    scalar_mult_su3_matrix_f(&Gamma, alpha2, &tmat);
    su3_adjoint_f(&tmat, lambda[dir] + i);

    /* the derivative which contributes to the globaln to the new global Sigma
       If this is the first level, then Sigma has to be initiallized. On later
       levels, we accumulate the respecive contributions
       */
    if (sigfresh == 0)
      scalar_mult_su3_matrix_f(&Gamma, alpha1, Sigma[dir] + i);
    else
      scalar_mult_sum_su3_matrix_f(&Gamma, alpha1, Sigma[dir] + i);
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute next-level Sigma (sig)
// lnk is the links the routine operates on
// dir1 is the direction of the link and the 'main' direction of Sigma
// lambda = Gamma^dag
void compute_sigma23(su3_matrix_f *sig, su3_matrix_f *lnk1, su3_matrix_f *lnk2,
                     su3_matrix_f *lambda1, su3_matrix_f *lambda2,
                     int dir1, int dir2) {

  register int i;
  register site *st;
  msg_tag *tag0, *tag1, *tag2, *tag3, *tag4;
  su3_matrix_f tmat, tmat2;

  // There are six terms: two staples each times three links

  /* get link[dir2] from direction dir1   */
  tag0 = start_gather_field(lnk2, sizeof(su3_matrix_f), dir1,
      EVENANDODD, gen_pt[0]);

  /* get Lambda[dir2] from direction dir1 */
  tag2 = start_gather_field(lambda2, sizeof(su3_matrix_f), dir1,
      EVENANDODD, gen_pt[2]);

  wait_gather(tag0);
  wait_gather(tag2);

  /* get link[dir1] from direction dir2   */
  tag1 = start_gather_field(lnk1, sizeof(su3_matrix_f), dir2,
      EVENANDODD, gen_pt[1]);

  /* get Lambda[dir1] from direction dir2 */
  tag3 = start_gather_field(lambda1, sizeof(su3_matrix_f), dir2,
      EVENANDODD, gen_pt[3]);

  // Prepare lower staple at x - dir2, store it in tempmat_nhyp2
  // and gather it to x
  FORALLSITES(i, st) {
    /* "term2" */
    mult_su3_nn_f(lambda1+i, (su3_matrix_f *)gen_pt[0][i], &tmat);
    mult_su3_an_f(&tmat, lnk2+i, tempmat_nhyp2+i);
    /* "term3" */
    mult_su3_nn_f(lnk1+i, (su3_matrix_f *)gen_pt[2][i], &tmat);
    mult_su3_an_f(&tmat, lnk2+i, &tmat2);
    add_su3_matrix_f(tempmat_nhyp2 + i, &tmat2, tempmat_nhyp2 + i);
    /* "term4" */
    mult_su3_nn_f(lnk1 + i, (su3_matrix_f *)gen_pt[0][i], &tmat);
    mult_su3_an_f(&tmat, lambda2 + i, &tmat2);
    add_su3_matrix_f(tempmat_nhyp2 + i, &tmat2, tempmat_nhyp2 + i);
  }

  /* gather staple from direction -dir2 to "home" site */
  tag4 = start_gather_field(tempmat_nhyp2, sizeof(su3_matrix_f),
      OPP_DIR(dir2), EVENANDODD, gen_pt[4]);

  wait_gather(tag1);
  wait_gather(tag3);

  /* Upper staple */
  FORALLSITES(i, st) {
    /* "term1" */
    mult_su3_nn_f(lambda2+i, (su3_matrix_f *)gen_pt[1][i], &tmat);
    mult_su3_na_f((su3_matrix_f *)gen_pt[0][i], &tmat, sig+i);
    /* "term5" */
    mult_su3_na_f((su3_matrix_f *)gen_pt[2][i],
        (su3_matrix_f *)gen_pt[1][i], &tmat);
    mult_su3_na_f(&tmat, lnk2+i, &tmat2);
    add_su3_matrix_f(sig+i, &tmat2, sig+i);
    /* "term6" */
    mult_su3_na_f((su3_matrix_f *)gen_pt[0][i],
        (su3_matrix_f *)gen_pt[3][i], &tmat);
    mult_su3_na_f(&tmat, lnk2+i, &tmat2);
    add_su3_matrix_f(sig+i, &tmat2, sig+i);
  }

  /* finally add the lower staple. */
  wait_gather(tag4);

  FORALLSITES(i, st)
    add_su3_matrix_f(sig+i, (su3_matrix_f *)gen_pt[4][i], sig+i);

  cleanup_gather(tag0);
  cleanup_gather(tag1);
  cleanup_gather(tag2);
  cleanup_gather(tag3);
  cleanup_gather(tag4);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Third-level force
#if  SMEAR_LEVEL == 3
void nhyp_force3(int dir3, int dir2) {
  register int i, dir, dir1, dir4;
  register site *s;

  FORALLUPDIR(dir) {
    if (dir == dir2 || dir == dir3)
      continue;
    FORALLUPDIR(dir4) {
      if (dir4 != dir && dir4 != dir3 && dir4 != dir2)
        break;
    }
    Sigma_update1(dir, SigmaH[dir], Staple1[dir4][dir], Lambda2,
                  1.0 - alpha_smear[2], alpha_smear[2] / 2.0, 1);
  }

  FORALLUPDIR(dir) {
    if (dir == dir2 || dir == dir3)
      continue;
    FORALLUPDIR(dir1) {
      if (dir1 == dir || dir1 == dir2 || dir1 == dir3)
        continue;
      compute_sigma23(tempmat_nhyp1, gauge_field_thin[dir],
                      gauge_field_thin[dir1], Lambda2[dir], Lambda2[dir1],
                      dir, dir1);
      FORALLSITES(i, s)
        sum_su3_matrix_f(tempmat_nhyp1 + i, Sigma[dir] + i);
    }
  }
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Second-level force
#if SMEAR_LEVEL > 1
void nhyp_force2(int dir2) {
  register int i;
  register site *s;
  int dir, dir3, dir4;
#if SMEAR_LEVEL == 3
  int imap[4][4] = {{0, 0, 1, 2}, {0, 0, 0, 3}, {1, 0, 0, 0}, {2, 3, 0, 0}};
  int iimap;
#endif

  // dir3 is the main direction of the twice-smeared hyplink2
  // dir2 is the secondary direction of the twice-smeared hyplink2
  // dir is the main direction of the once-smeared link
  FORALLUPDIR(dir3) {
    if (dir3 != dir2) {
      Sigma_update1(dir3, SigmaH[dir3], Staple2[dir2][dir3], Lambda1,
                    1.0 - alpha_smear[1], alpha_smear[1] / 4.0, 1);
    }
  }

  FORALLUPDIR(dir3) {
    if (dir3 == dir2)
      continue;
    FORALLUPDIR(dir) {
      if (dir3 == dir || dir == dir2)
        continue;
      FORALLUPDIR(dir4) {
        if (dir4 != dir && dir4 != dir2 && dir4 != dir3)
          break;
      }
#if SMEAR_LEVEL == 3
      compute_sigma23(SigmaH[dir], hyplink1[dir4][dir], hyplink1[dir4][dir3],
                      Lambda1[dir], Lambda1[dir3], dir, dir3);
#else   // SMEAR_LEVEL == 2
      compute_sigma23(SigmaH[dir], gauge_field_thin[dir],
                      gauge_field_thin[dir3], Lambda1[dir], Lambda1[dir3],
                      dir, dir3);

      FORALLSITES(i, s)
        sum_su3_matrix_f(SigmaH[dir] + i, Sigma[dir] + i);
#endif
    }

    // This part is really awkward
    // nhyp_force3 is symmetric in the arguments, but the input is not
    // So add the dir2, dir3 and the dir3, dir2 terms
    // To do so, one has to sort them...
    // To save memory, store only the upper triangular part of the 4x4 matrix
    // For that, only a 4 'vector' is necessary
    // Which field to use is given by the imap array hard-coded above
#if SMEAR_LEVEL == 3
    iimap = imap[dir2][dir3];
    if (dir2 < dir3) {
      FORALLSITES(i, s) {
        FORALLUPDIR(dir)
          su3mat_copy_f(SigmaH[dir] + i, SigmaH2[iimap][dir] + i);
      }
    }
    else {
      FORALLSITES(i, s) {
        FORALLUPDIR(dir)
          sum_su3_matrix_f(SigmaH2[iimap][dir] + i, SigmaH[dir] + i);
      }
      nhyp_force3(dir2, dir3);
    }
#endif
  }
}
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// First-level force
void nhyp_force1() {
  int dir, dir2;
#if SMEAR_LEVEL == 1
  register int i;
  register site *s;
#endif

  // Loop over the link directions,
  // compute sigma from the link itself,
  // construct the Lambdas
  FORALLUPDIR(dir) {
    Sigma_update1(dir, Sigma[dir], Staple3[dir], LambdaU,
                  1.0 - alpha_smear[0], alpha_smear[0] / 6.0, 0);
  }

  // Construct field_offsets pointing to the links in dir2 direction
  // where dir1 is excluded.  Here, this makes no sense, however this
  // mechanism makes compute_sigma23 re-useable on each level
  FORALLUPDIR(dir2) {
    FORALLUPDIR(dir) {
      if (dir == dir2)
        continue;
#if SMEAR_LEVEL > 1
      compute_sigma23(SigmaH[dir], hyplink2[dir2][dir], hyplink2[dir][dir2],
                      LambdaU[dir], LambdaU[dir2], dir, dir2);
#else  // SMEAR_LEVEL == 1
      compute_sigma23(SigmaH[dir], gauge_field_thin[dir],
                      gauge_field_thin[dir2], LambdaU[dir], LambdaU[dir2],
                      dir, dir2);

      FORALLSITES(i, s)
        sum_su3_matrix_f(SigmaH[dir] + i, Sigma[dir] + i);
#endif
    }
#if SMEAR_LEVEL > 1
    nhyp_force2(dir2);
#endif
  }
}
// -----------------------------------------------------------------
