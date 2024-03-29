// -----------------------------------------------------------------
// Fix Coulomb or Lorenz gauge by doing successive SU(2) gauge hits
// Double-precision global sums, as always
// Reunitarize after a pre-set interval defined below

// Prototype:
// void gaugefix(int gauge_dir, Real relax_boost, int max_gauge_iter,
//               Real gfix_tol, field_offset diffmat, field_offset sumvec);
//
// EXAMPLE: Fixing only the fundamental link matrices to Coulomb gauge,
//          using scratch space in the site structure
//          (matrix pm and vector chi)
// gaugefix(TUP, 1.5, 500, 1.0e-7, F_OFFSET(mp), F_OFFSET(chi));
//
// gauge_dir is the "time" direction used to define Coulomb gauge:
//   TUP for evaluating propagators in the time-like direction
//   ZUP for screening lengths
//   8 (or anything except 0, 1, 2, 3) for Lorenz gauge
//   (cf. definition of FORALLUPDIRBUT in ../include/macros.h)
// relax_boost is the overrelaxation parameter
// max_gauge_iter is the maximum number of gauge-fixing iterations
// gfix_tol tells us to stop if the action change is less than this
// diffmat is scratch space for a matrix
// sumvec is scratch space for a vector
// If diffmat or sumvec are negative, scratch space is malloced here
#include "generic_includes.h"
#define REUNIT_INTERVAL 20

// Scratch space
static matrix *diffmatp;               // Malloced diffmat pointer
static vector *sumvecp;                // Malloced sumvec pointer
field_offset diffmat_offset, sumvec_offset;   // Field offsets
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Accumulate sums and differences of link matrices
// for determining optimum hit for gauge fixing
// Differences are kept in diffmat
// Diagonal elements of the sums are kept in sumvec
// Downward links gathered to gen_pt[dir] in gaugefixstep()
void accum_gauge_hit(int gauge_dir, int parity) {
  register int i, j, dir;
  register matrix *m1, *m2;
  register site *s;

  // Clear sumvec and diffmat
  FORSOMEPARITY(i, s, parity) {
    if (diffmat_offset >= 0)
      clear_mat((matrix *)F_PT(s, diffmat_offset));
    else
      clear_mat(&diffmatp[i]);
    if (sumvec_offset >= 0)
      clearvec((vector *)F_PT(s, sumvec_offset));
    else
      clearvec(&sumvecp[i]);
  }

  // Subtract upward link contributions
  FORSOMEPARITY(i, s, parity) {
    FORALLUPDIRBUT(gauge_dir, dir) {
      // Upward link matrix
      m1 = &(s->link[dir]);
      if (diffmat_offset >= 0)
        dif_mat(m1, (matrix *)F_PT(s, diffmat_offset));
      else
        dif_mat(m1, &diffmatp[i]);

      if (sumvec_offset >= 0) {
        for (j = 0; j < NCOL; j++)
          CSUM(((vector *)F_PT(s, sumvec_offset))->c[j], m1->e[j][j]);
      }
      else {
        for (j = 0; j < NCOL; j++)
          CSUM(sumvecp[i].c[j], m1->e[j][j]);
      }
    }
  }

  // Add downward link contributions
  FORSOMEPARITY(i, s, parity) {
    FORALLUPDIRBUT(gauge_dir, dir) {
      // Downward link matrix -- gathered in gaugefixstep()
      m2 = (matrix *)gen_pt[dir][i];

      if (diffmat_offset >= 0)
        sum_mat(m2, (matrix *)F_PT(s, diffmat_offset));
      else
        sum_mat(m2, &diffmatp[i]);

      if (sumvec_offset >= 0) {
        for (j = 0; j < NCOL; j++)
          CSUM(((vector *)F_PT(s, sumvec_offset))->c[j], m2->e[j][j]);
      }
      else {
        for (j = 0; j < NCOL; j++)
          CSUM(sumvecp[i].c[j], m2->e[j][j]);
      }
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Do optimum SU(2) gauge hit for p, q subspace
void do_hit(int gauge_dir, int parity, int p, int q, Real relax_boost) {
  register int dir, i;
  register site *s;
  Real a0, a1, a2, a3, asq, a0sq, x, r, xdr;
  su2_matrix u;

  // Accumulate sums for determining optimum gauge hit
  accum_gauge_hit(gauge_dir, parity);

  FORSOMEPARITY(i, s, parity) {
    // The SU(2) hit matrix is represented as a0 + i * Sum j (sigma j * aj)
    // The locally optimum unnormalized components a0, aj are determined
    // from the current link in direction dir and the link downlink
    // in the same direction on the neighbor in the direction opposite dir
    // The expression is
    //   a0 = Sum dir Tr Re 1       * (downlink dir + link dir)
    //   aj = Sum dir Tr Im sigma j * (downlink dir - link dir), j = 1, 2, 3
    // where 1, sigma j are unit and Pauli matrices on the p, q subspace
    // a0 =  s->sumvec.c[p].real + s->sumvec.c[q].real;
    // a1 =  s->diffmat.e[q][p].imag + s->diffmat.e[p][q].imag;
    // a2 = -s->diffmat.e[q][p].real + s->diffmat.e[p][q].real;
    // a3 =  s->diffmat.e[p][p].imag - s->diffmat.e[q][q].imag;
    if (sumvec_offset >= 0) {
      a0 = ((vector *)F_PT(s, sumvec_offset))->c[p].real
         + ((vector *)F_PT(s, sumvec_offset))->c[q].real;
    }
    else
      a0 = sumvecp[i].c[p].real + sumvecp[i].c[q].real;

    if (diffmat_offset >= 0) {
      a1 =  ((matrix *)F_PT(s, diffmat_offset))->e[q][p].imag
         +  ((matrix *)F_PT(s, diffmat_offset))->e[p][q].imag;
      a2 = -((matrix *)F_PT(s, diffmat_offset))->e[q][p].real
         +  ((matrix *)F_PT(s, diffmat_offset))->e[p][q].real;
      a3 =  ((matrix *)F_PT(s, diffmat_offset))->e[p][p].imag
         -  ((matrix *)F_PT(s, diffmat_offset))->e[q][q].imag;
    }
    else {
      a1 =  diffmatp[i].e[q][p].imag + diffmatp[i].e[p][q].imag;
      a2 = -diffmatp[i].e[q][p].real + diffmatp[i].e[p][q].real;
      a3 =  diffmatp[i].e[p][p].imag - diffmatp[i].e[q][q].imag;
    }

    // Over-relaxation boost
    // This algorithm is designed to give little change for large |a|
    // and to scale up the gauge transformation by a factor of relax_boost
    // for small |a|
    asq = a1 * a1 + a2 * a2 + a3 * a3;
    a0sq = a0 * a0;
    x = (relax_boost * a0sq + asq) / (a0sq + asq);
    r = sqrt((double)(a0sq + x * x * asq));
    xdr = x / r;

    // Normalize and boost
    a0 = a0 / r;
    a1 = a1 * xdr;
    a2 = a2 * xdr;
    a3 = a3 * xdr;

    // Elements of SU(2) matrix
    u.e[0][0] = cmplx(a0, a3);
    u.e[0][1] = cmplx(a2, a1);
    u.e[1][0] = cmplx(-a2, a1);
    u.e[1][1] = cmplx(a0, -a3);

    // Do SU(2) hit on all upward links
    FORALLUPDIR(dir)
      left_su2_hit_n(&u, p, q,&(s->link[dir]));

    // Do SU(2) hit on all downward links
    FORALLUPDIR(dir)
      right_su2_hit_a(&u, p, q,(matrix *)gen_pt[dir][i]);
  }
  // Exit with modified downward links left in communications buffer
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Return average of gauge fixing action for sites of the given parity
// Normalize to a maximum of 1 when all links are unit matrices
double get_gauge_fix_action(int gauge_dir, int parity) {
  register int dir, i, ndir = 0;
  register site *s;
  double gauge_fix_action = 0.0;
  complex tc;

  FORSOMEPARITY(i, s, parity) {
    FORALLUPDIRBUT(gauge_dir, dir) {
      tc = trace(&(s->link[dir]));
      gauge_fix_action += (double)tc.real;

      tc = trace((matrix *)gen_pt[dir][i]);
      gauge_fix_action += (double)tc.real;
    }
  }

  // Count number of terms to average
  FORALLUPDIRBUT(gauge_dir, dir)
    ndir++;

  // Sum over all sites of this parity
  g_doublesum(&gauge_fix_action);

  // Average is normalized to max of 1 / 2 on sites of one parity
  gauge_fix_action *= one_ov_N / ((double)(2.0 * ndir * volume));
  return gauge_fix_action;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Carries out one iteration in the gauge-fixing process
void gaugefixstep(int gauge_dir, double *av_gauge_fix_action,
                  Real relax_boost) {

  register int dir, i, c1, c2;
  register site *s;
  int parity;
  msg_tag *mtag[8];
  Real gauge_fix_action;

  // Alternate parity to prevent interactions during gauge transformation
  *av_gauge_fix_action = 0.0;
  g_sync();
  fflush(stdout);

  for (parity = ODD; parity <= EVEN; parity++) {
    // Gather all downward links at once
    FORALLUPDIR(dir) {
      mtag[dir] = start_gather_site(F_OFFSET(link[dir]), sizeof(matrix),
                                    OPP_DIR(dir), parity, gen_pt[dir]);
    }
    FORALLUPDIR(dir)
      wait_gather(mtag[dir]);

    // Do optimum gauge hit on all subspaces
    for (c1 = 0; c1 < NCOL - 1; c1++) {
      for (c2 = c1 + 1; c2 < NCOL; c2++)
        do_hit(gauge_dir, parity, c1, c2, relax_boost);
    }

    // Total gauge fixing action for sites of this parity
    gauge_fix_action = get_gauge_fix_action(gauge_dir, parity);
    *av_gauge_fix_action += gauge_fix_action;

    // Scatter downward link matrices by gathering to sites of opposite parity
    FORALLUPDIR(dir) {
      // Synchronize before scattering to be sure the new modified link
      // matrices are all ready to be scattered and diffmat is not
      // overwritten before it is used
      g_sync();

      // First copy modified link for this dir
      // from comm buffer or node to diffmat
      FORSOMEPARITY(i, s, parity) {
        if (diffmat_offset >= 0) {
          mat_copy((matrix *)(gen_pt[dir][i]),
                   (matrix *)F_PT(s, diffmat_offset));
        }
        else
          mat_copy((matrix *)(gen_pt[dir][i]), &diffmatp[i]);
      }
      cleanup_gather(mtag[dir]);    // Done with gen_pt[dir]

      // Synchronize to make sure the previous copy happens before the
      // subsequent gather below
      g_sync();

      // Gather diffmat onto sites of opposite parity
      if (diffmat_offset >= 0) {
        mtag[dir] = start_gather_site(diffmat_offset, sizeof(matrix),
                                      dir, OPP_PAR(parity), gen_pt[dir]);
      }
      else {
        mtag[dir] = start_gather_field(diffmatp, sizeof(matrix),
                                       dir, OPP_PAR(parity), gen_pt[dir]);
      }
      wait_gather(mtag[dir]);

      // Copy modified matrices into proper location
      FORSOMEPARITY(i, s, OPP_PAR(parity))
        mat_copy((matrix *)(gen_pt[dir][i]), &(s->link[dir]));

      cleanup_gather(mtag[dir]);
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void gaugefixscratch(field_offset diffmat, field_offset sumvec) {
  diffmat_offset = diffmat;
  diffmatp = NULL;
  if (diffmat_offset < 0) {
    diffmatp = malloc(sites_on_node * sizeof(*diffmatp));
    if (diffmatp == NULL) {
      node0_printf("gaugefix: Can't malloc diffmat\n");
      fflush(stdout);
      terminate(1);
    }
  }

  sumvec_offset = sumvec;
  sumvecp = NULL;
  if (sumvec_offset < 0) {
    sumvecp = malloc(sites_on_node * sizeof(*sumvecp));
    if (sumvecp == NULL) {
      node0_printf("gaugefix: Can't malloc sumvec\n");
      fflush(stdout);
      terminate(1);
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void gaugefix(int gauge_dir, Real relax_boost, int max_gauge_iter,
              Real gfix_tol, field_offset diffmat, field_offset sumvec) {

  int gauge_iter;
  double current_av, old_av = 0.0, del_av = 0.0;

  // Require at least 8 gen_pt values for gauge fixing
  if (N_POINTERS < 8) {
    node0_printf("gaugefix: N_POINTERS must be at least 8\n");
    fflush(stdout);
    terminate(1);
  }

  // Set up work space
  gaugefixscratch(diffmat, sumvec);

  // Do at most max_gauge_iter iterations
  // Stop after second step if the change in the average gauge fixing action
  // is smaller than gfix_tol
  for (gauge_iter = 0; gauge_iter < max_gauge_iter; gauge_iter++) {
    gaugefixstep(gauge_dir, &current_av, relax_boost);

    if (gauge_iter != 0) {
      del_av = current_av - old_av;
      if (fabs(del_av) < gfix_tol)
        break;
    }
    old_av = current_av;

    // Reunitarize when iteration count is a multiple of REUNIT_INTERVAL
    if ((gauge_iter % REUNIT_INTERVAL) == (REUNIT_INTERVAL - 1)) {
      reunitarize();
#ifdef DEBUG_CHECK
      // Print upon reunitarization
      node0_printf("step %d av gf action %.6g, delta %.3g\n",
                   gauge_iter, current_av, del_av);
#endif
    }
  }
  // Reunitarize at the end, unless we just did it in the loop
  if ((gauge_iter % REUNIT_INTERVAL) != 0)
    reunitarize();

  // Free workspace
  if (diffmat_offset < 0)
    free(diffmatp);
  if (sumvec_offset < 0)
    free(sumvecp);

  node0_printf("GFIX ended at step %d. Ave gf action %.6g, delta %.3g\n",
               gauge_iter, current_av, del_av);
  if (gauge_iter == max_gauge_iter)
    node0_printf("WARNING: max_gauge_iter %d saturated\n", gauge_iter);
}
// -----------------------------------------------------------------
