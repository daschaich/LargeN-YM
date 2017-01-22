// -----------------------------------------------------------------
// Solve (Mdag.M) psi = chi for clover fermions in dynamical HMC updates

// Use the LU preconditioned fermion matrix
// Solution and result live only on one parity of the lattice
// If dslash_oe takes a source on even sites to a result on odd sites, then
//   M = R_e - kappa^2 * dslash_eo * R_o^(-1) * dslash_oe

// Restart after each niter passes, up to CG_RESTART times
#define CG_RESTART 0

// rsqmin = desired rsq, quit when we reach rsq = rsqmin*source_norm.
#include "cl_dyn_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Wrapper for (unshifted) fermion operator and its adjoint
// CG uses special restarted dslashes, which are more efficient than this
// The sign is PLUS for M, MINUS for MDAG
// !!! Assume make_clov(CKU0) and make_clovinv(ODD) have been called
//   The first computes R = 1 - i*CKU0 sigma_{\mu\nu} F_{\mu\nu} on each site
//   The second inverts R_o on odd sites only
//   This fixes the operator to the EVEN sublattice
//   (So we can use *src = *dest if we want)
// Uses tempwvec and p for temporary storage, preserving src
void fermion_op(wilson_vector *src, wilson_vector *dest,
                int sign, int parity) {

  register int i, otherparity = ODD;
  register site *s;

  switch(parity) {
    case EVEN:       otherparity = ODD;        break;
    case ODD:        otherparity = EVEN;       break;
    default:
      node0_printf("fermion_op parity can be only EVEN or ODD\n");
      terminate(1);
  }

  // M = R_p - kappa^2 * dslash_po * R_o^(-1) * dslash_op
  // 1) tempwvec_e <-- R_e src_e, then untouched until last step
  mult_ldu(src, tempwvec, parity);

  // 2) tempwvec_o <-- dslash_oe src_e
  dslash(src, tempwvec, sign, otherparity);

  // 3) p_o <-- R_o^(-1) tempwvec_o
  mult_ldu(tempwvec, p, otherparity);

  // 4) dest_e <-- dslash_eo p_o
  dslash(p, dest, sign, parity);

  // 5) dest_e <-- tempwvec_e - kappa^2 * chi[0]_e and done
  FORSOMEPARITY(i, s, parity) {
    scalar_mult_wvec(&(dest[i]), mkappaSq, &(dest[i]));
    sum_wvec(&(tempwvec[i]), &(dest[i]));
  }
}

// As above, but with restarted gathers where possible
void fermion_op_special(wilson_vector *src, wilson_vector *dest,
                        int sign, int parity, msg_tag **o_tag,
                        msg_tag **p_tag, int is_started) {

  register int i, otherparity = ODD;
  register site *s;

  switch(parity) {
    case EVEN:       otherparity = ODD;        break;
    case ODD:        otherparity = EVEN;       break;
    default:
      node0_printf("fermion_op parity can be only EVEN or ODD\n");
      terminate(1);
  }

  mult_ldu(src, tempwvec, parity);
  dslash_special(src, tempwvec, sign, otherparity, o_tag, is_started);
  mult_ldu(tempwvec, p, otherparity);
  dslash_special(p, dest, sign, parity, p_tag, is_started);
  FORSOMEPARITY(i, s, parity) {
    scalar_mult_wvec(&(dest[i]), mkappaSq, &(dest[i]));
    sum_wvec(&(tempwvec[i]), &(dest[i]));
  }
}
// -----------------------------------------------------------------




// -----------------------------------------------------------------
// The level tells us which source and initial guess/answer to use:
//   chi[level] and psi[level], respectively
// mshift is the Hasenbusch shift
// Uses temporaries r (residual), p, mp, and tempwvec
// Return total number of Mdag.M applications (summed over restarts)
int congrad(int level, Real mshift, int parity) {
  register int i;
  register site *s;
  int iter = 0;
  Real a, b, MSq = mshift * mshift;
  double rsqstop = 0.0, rsq = 0.0, oldrsq = 0.0;
  double source_norm = 0.0, pkp = 0.0;      // pkp = p Mdag.M p
  msg_tag *tag[8], *tag2[8];
#ifdef DEBUG_CHECK
  double mflops, dtime = -dclock();
#endif

#ifdef TIMING
  TIC(0)
#endif

#ifdef DEBUG_CHECK
  printf("Dumping source...\n");
  FORSOMEPARITY(i, s, parity) {
    printf("chi[%d](%d, %d, %d, %d):\n", level, s->x, s->y, s->z, s->t);
    dump_wvec(&(chi[level][i]));
  }
#endif

start:
  // 1) mp <-- Mdag.M on psi[level]
  // 2) r, p <-- chi[level] - mp
  // 3) rsq = |r|^2
  // 4) source_norm = |chi[level]|^2
  rsq = 0.0;
  source_norm = 0.0;

  // 1) mp <-- Mdag.M on psi[level]
  fermion_op_special(psi[level], mp, PLUS, parity, tag, tag2, 0);
  fermion_op_special(mp, mp, MINUS, parity, tag, tag2, 1);
  FORSOMEPARITY(i, s, parity) {
    // Add mshift^2 psi[level] to complete Mdag.M on psi[level]
    if (mshift > 0.0)
      scalar_mult_sum_wvec(&(psi[level][i]), MSq, &(mp[i]));

    // 2) r, p <-- chi[level] - mp
    sub_wvec(&(chi[level][i]), &(mp[i]), &(r[i]));
    copy_wvec(&(r[i]), &(p[i]));

    // 3) rsq = |r|^2
    magsq_wvec_sum(&(r[i]), &rsq);

    // 4) source_norm = |chi[level]|^2
    magsq_wvec_sum(&(chi[level][i]), &source_norm);
  }
  g_doublesum(&source_norm);
  g_doublesum(&rsq);
  iter++;                       // Count applications of Mdag.M

  rsqstop = rsqmin * source_norm;
#ifdef CG_DEBUG
  node0_printf("CG source_norm = %.4g ", source_norm);
  node0_printf("--> rsqstop = %.4g\n", rsqstop);
  node0_printf("CG iter %d, rsq %.4g, pkp %.4g, a %.4g\n",
               iteration, rsq, pkp, a);
#endif
  if (rsq <= rsqstop) {
    total_iters += iter;        // Global count of Mdag.M applications
    FORALLUPDIR(i) {
      cleanup_gather(tag[i]);
      cleanup_gather(tag2[i]);
      cleanup_gather(tag[OPP_DIR(i)]);
      cleanup_gather(tag2[OPP_DIR(i)]);
    }

#ifdef TIMING
    TOC(0, time_dcongrad)
#endif

    return iter;
  }

  // Main loop -- do until convergence or time to restart
  // 1) oldrsq <-- rsq
  // 2) mp <-- Mdag.M p on given parity
  // 3) pkp <-- sum(p Mdag.M p)
  // 4) a <-- rsq / pkp
  // 5) psi[level] <-- psi[level] + a * p
  // 6) r <-- r - a * mp
  // 7) rsq <-- sum(|r|^2)
  // 8) b <-- rsq / oldrsq
  // 9) p <-- r + b * p
  do {    // We accumulate iter potentially over several loops
    // 1) oldrsq <-- rsq
    oldrsq = rsq;

    // 2) mp <-- Mdag.M p on given parity
    fermion_op_special(p, mp, PLUS, parity, tag, tag2, 1);
    fermion_op_special(mp, mp, MINUS, parity, tag, tag2, 1);
    pkp = 0.0;
    FORSOMEPARITY(i, s, parity) {
      // Add mshift^2 p
      if (mshift > 0.0)
        scalar_mult_sum_wvec(&(p[i]), MSq, &(mp[i]));

      // 3) pkp <-- sum(p Mdag.M p)
      wvec_rdot_sum(&(p[i]), &(mp[i]), &pkp);
    }
    g_doublesum(&pkp);
    a = (Real)(rsq / pkp);      // 4) a <-- rsq / pkp
    iter++;                     // Count applications of Mdag.M

    rsq = 0.0;
    FORSOMEPARITY(i, s, parity) {
      // 5) psi[level] <-- psi[level] + a * p
      scalar_mult_sum_wvec(&(p[i]), a, &(psi[level][i]));

      // 6) r <-- r - a * mp
      scalar_mult_dif_wvec(&(mp[i]), a, &(r[i]));

      // 7) rsq <-- sum(|r|^2)
      magsq_wvec_sum(&(r[i]), &rsq);
    }
    g_doublesum(&rsq);
#ifdef CG_DEBUG
    node0_printf("CG iter %d, rsq %.4g, pkp %.4g, a %.4g\n",
                 iter, rsq, pkp, a);
#endif
    // See if we're done
    if (rsq <= rsqstop) {
      total_iters += iter;    // Global count of Mdag.M applications
      FORALLUPDIR(i) {
        cleanup_gather(tag[i]);
        cleanup_gather(tag2[i]);
        cleanup_gather(tag[OPP_DIR(i)]);
        cleanup_gather(tag2[OPP_DIR(i)]);
      }
#ifdef CG_DEBUG
      dtime += dclock();
      mflops = (double)(2840.0 * volume * iter / (1e6 * dtime * numnodes()));
      node0_printf("CG time = %.4g iters = %d mflops = %.4g\n",
                   dtime, iter, mflops);
#endif

#ifdef TIMING
      TOC(0, time_dcongrad)
#endif

      return iter;
    }

    // We're not done, so set up next iteration
    // 8) b <-- rsq / oldrsq
    b = rsq / oldrsq;
    FORSOMEPARITY(i, s, parity) {
      // 9) p <-- r + b * p
      scalar_mult_wvec(&(p[i]), b, &(p[i]));
      sum_wvec(&(r[i]), &(p[i]));
    }
  } while (iter % niter != 0);

  // We have done niter iterations, time to restart or give up
  total_iters += iter;    // Global count of Mdag.M applications
  FORALLUPDIR(i) {
    cleanup_gather(tag[i]);
    cleanup_gather(tag2[i]);
    cleanup_gather(tag[OPP_DIR(i)]);
    cleanup_gather(tag2[OPP_DIR(i)]);
  }

  // Hard-coded number of restarts in defines.h
  if (iter < CG_RESTART * niter)
    goto start;
  if (rsq > rsqstop) {
    node0_printf("WARNING: CG did not converge, size_r = %.2g\n",
                 sqrt(rsq / source_norm));
  }

#ifdef TIMING
  TOC(0, time_dcongrad)
#endif

  return iter;
}
// -----------------------------------------------------------------
