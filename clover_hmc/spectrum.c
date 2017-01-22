// -----------------------------------------------------------------
// Hadron spectrum for Wilson-clover fermions
// Seems to have NCOL = DIMF = 3 hard coded
#include "cl_dyn_includes.h"
#define POINT 1
#define WALL 2

#if (NCOL != 3 || DIMF != 3)
  #error "Spectrum only works if NCOL = DIMF = 3"
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void source(wilson_vector *src, int color, int spin, int type,
            int x0, int y0, int z0, int t0, Real gamma) {

  register int i;
  register site *s;
  short my_x, my_y, my_z;
  Real rx, ry, rz, radius2;

#ifdef DEBUG_CHECK
  node0_printf("source is making source type %d\n", type);
#endif

  // Clear source
  FORALLSITES(i, s)
    clear_wvec(&(src[i]));

  if (type == POINT) {
    // Unit delta function at (x0, y0, z0, t0)
    if (node_number(x0, y0, z0, t0) == mynode()) {
      i = node_index(x0, y0, z0, t0);
      src[i].d[spin].c[color].real = 1.0;
    }
  }
  else if (type == WALL) {
    // Gaussian fixed to timeslice t0, centered on (x0, y0, z0)
    FORALLSITES(i, s) {
      if (s->t != t0)
        continue;

      my_x = ((s->x) - x0 + nx) % nx;
      rx = (my_x < (nx - my_x)) ? (Real)my_x : (Real)(my_x - nx);
      my_y = ((s->y) - y0 + ny) % ny;
      ry = (my_y < (ny - my_y)) ? (Real)my_y : (Real)(my_y - ny);
      my_z = ((s->z) - z0 + nz) % nz;
      rz = (my_z < (nz - my_z)) ? (Real)my_z : (Real)(my_z - nz);

      radius2 = rx * rx + ry * ry + rz * rz;
      src[i].d[spin].c[color].real = (Real)exp((double)(-radius2 * gamma));
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Pion and rho with point sink and zero spatial momentum
// Sources are type wilson_matrix
void meson(field_offset src1, field_offset src2, char id_string[]) {
  register int i, t;
  register site *s;
  int my_t, cf, sf, ci, si;
  // These are time-slice sums of propagators,
  // with and without an extra gamma0 at the sink
  Real *prop = malloc(nt * sizeof(*prop));
  Real *prop0 = malloc(nt * sizeof(*prop0));
  Real *propim = malloc(nt * sizeof(*propim));
  Real *prop0im = malloc(nt * sizeof(*prop0im));
  complex g1, g2;
  wilson_matrix tmat, tmat2;       // Temporary storage

  // Measure the pion propagators -- gamma5 gamma5 correlator and
  //                                 gamma5 gamma5gamma0 correlator
  for (t = 0; t < nt; t++) { /* clear meson propgators */
    prop[t] = 0.0;
    prop0[t] = 0.0;
    propim[t] = 0.0;
    prop0im[t] = 0.0;
  }

  /*first, dirac multiplication by the source gamma matrices */
  FORALLSITES(i, s) {
    my_t = s->t;
    /*first, dirac multiplication by the source gamma matrices (on left) */
    mult_by_gamma_left(((wilson_matrix *)F_PT(s, src1)), &tmat2,
                       GAMMAFIVE);
    /*first, dirac multiplication by gamma-5 (antiquark will need this
      along with complex conjugation */
    mult_by_gamma_left(&tmat2, &tmat, GAMMAFIVE);

    /* dirac multiplication by the sink gamma matrices (on right) */
    /** Commented out for this example, because gamma5^2 = 1 **
      mult_by_gamma_right(&tmat, &tmat2, GAMMAFIVE);
      mult_by_gamma_right(&tmat2, &tmat, GAMMAFIVE);
     **/
    /* Multiply by gamma0 for gamma5 - gamma5*gamma0 correlator */
    mult_by_gamma_right(&tmat, &tmat2, TUP);

    // Trace over propagators
    for (si = 0; si < 4; si++) {
      for (ci = 0; ci < 3; ci++) {
        for (sf = 0; sf < 4; sf++) {
          for (cf = 0; cf < 3; cf++) {
            g2 = ((wilson_matrix *)F_PT(s, src2))->d[si].c[ci].d[sf].c[cf];
            g1 = tmat.d[si].c[ci].d[sf].c[cf];
            prop[my_t] += g1.real * g2.real + g1.imag * g2.imag;
            propim[my_t] += g1.imag * g2.real - g1.real * g2.imag;
            g1 = tmat2.d[si].c[ci].d[sf].c[cf];
            prop0[my_t] += g1.real * g2.real + g1.imag * g2.imag;
            prop0im[my_t] += g1.imag * g2.real - g1.real * g2.imag;
          }
        }
      }
    }
  }
  // Print the pion propagators
  for (t = 0; t < nt; t++) {
    g_floatsum(prop + t);
    g_floatsum(prop0 + t);
    g_floatsum(propim + t);
    g_floatsum(prop0im + t);
    node0_printf("%sPSEUDO2 %d %.8g %.8g %.8g %.8g\n", id_string, t,
                 prop[t], propim[t], prop0[t], prop0im[t]);
  }

  // Measure the rho propagator -- gamma0gamma3 gamma0gamma3 correlator
  // Clear meson propagators
  for (t = 0; t < nt; t++) {
    prop[t] = 0.0;
    propim[t] = 0.0;
  }
  FORALLSITES(i, s) {
    my_t = s->t;
    /*first, dirac multiplication by the source gamma matrices (on left) */
    mult_by_gamma_left(((wilson_matrix *)F_PT(s, src1)), &tmat2, ZUP);
    mult_by_gamma_left(&tmat2, &tmat, TUP);
    /* dirac multiplication by gamma5 (antiquark will need this
       along with complex conjugation) */
    mult_by_gamma_left(&tmat, &tmat2, GAMMAFIVE);

    /* dirac multiplication by the sink gamma matrices (on right) */
    mult_by_gamma_right(&tmat2, &tmat, TUP);
    mult_by_gamma_right(&tmat, &tmat2, ZUP);
    /* dirac multiplication by gamma5 (finishing up antiquark) */
    mult_by_gamma_right(&tmat2, &tmat, GAMMAFIVE);

    // Trace over propagators
    for (si = 0; si < 4; si++) {
      for (ci = 0; ci < 3; ci++) {
        for (sf = 0; sf < 4; sf++) {
          for (cf = 0; cf < 3; cf++) {
            g1 = tmat.d[si].c[ci].d[sf].c[cf];
            g2 = ((wilson_matrix *)F_PT(s, src2))->d[si].c[ci].d[sf].c[cf];
            /*minus sign since we should be using adjoint of
               hadron operator at one end */
            prop[my_t] -= g1.real * g2.real + g1.imag * g2.imag;
            propim[my_t] -= g1.imag * g2.real - g1.real * g2.imag;
          }
        }
      }
    }
  }
  // Print the rho propagator
  for (t = 0; t < nt; t++) {
    g_floatsum(prop + t);
    g_floatsum(propim + t);
    node0_printf("%sRHO03032 %d %.8g %.8g\n",
                 id_string, t, (double)prop[t], (double)propim[t]);
  }

  // Measure the rho propagator -- gamma3 gamma3 correlator
  // Clear meson propagators
  for (t = 0; t < nt; t++) {
    prop[t] = 0.0;
    propim[t] = 0.0;
  }
  FORALLSITES(i, s) {
    my_t = s->t;
    /*first, dirac multiplication by the source gamma matrices (on left) */
    mult_by_gamma_left(((wilson_matrix *)F_PT(s, src1)), &tmat2, ZUP);
    /*first, dirac multiplication by gamma5 (antiquark will need this
      along with complex conjugation) */
    mult_by_gamma_left(&tmat2, &tmat, GAMMAFIVE);

    /* dirac multiplication by the sink gamma matrices (on right) */
    mult_by_gamma_right(&tmat, &tmat2, ZUP);
    /* dirac multiplication by gamma5 (finishing up antiquark) */
    mult_by_gamma_right(&tmat2, &tmat, GAMMAFIVE);

    // Trace over propagators
    for (si = 0; si < 4; si++) {
      for (ci = 0; ci < 3; ci++) {
        for (sf = 0; sf < 4; sf++) {
          for (cf = 0; cf < 3; cf++) {
            g1 = tmat.d[si].c[ci].d[sf].c[cf];
            g2 = ((wilson_matrix *)F_PT(s, src2))->d[si].c[ci].d[sf].c[cf];
            prop[my_t] += g1.real * g2.real + g1.imag * g2.imag;
            propim[my_t] += g1.imag * g2.real - g1.real * g2.imag;
          }
        }
      }
    }
  }
  // Print the rho propagator
  for (t = 0; t < nt; t++) {
    g_floatsum(prop + t);
    g_floatsum(propim + t);
    node0_printf("%sRHO332 %d %.8g %.8g\n", id_string, t,
                 (double)(*(prop + t)), (double)(*(propim + t)));
  }
  free(prop);
  free(prop0);
  free(propim);
  free(prop0im);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
/* Baryon propagator:
 propagator            point to plane baryon correlation (write)

 We are using a chiral basis for Dirac matrices (gamma5 diagonal)
 In that basis relativistic wave functions are much easier to implement.
 For the proton we have
         |p, ispinp> = (u C gamma5 d)u_ispinp     (1)
                     = (u_1 d_2 - u_2 d_1 + u_3 d_4 -u_4 d_3)u_ispinp
 for a proton of Dirac index (= helicity in a chiral basis).
  The wave function is encoded in the terms chi_in(i, j, k) and chi_out(i, j, k)
 where the first index labels the spinor of the d quark and the
 last two indices are for the two u quarks.
 The Delta is similar, only all plus signs

  We use the unsymmetrized w. f. and
  include direct and exchange terms in the propagator explicitly
         ___
         \
 b(t) =   >  < b(t , 0) b(t + t, x) >
         /        1        1
         ---
          x
*/

#define Nc 3
#define Ns 4

#define PROTON 0
#define DELTA 1
#define DELTA_N 2
// Sources are type wilson_matrix
// ispinp1 and ispinp2 are the spins of the source and sink -- always zero
void baryon(field_offset src1, field_offset src2, field_offset src3) {
  register int isite, t;
  register site *s;
  int my_t, i, j, k, ispinp1 = 0, ispinp2 = 0;
  int si_1, si_2, si_3, sf_1, sf_2, sf_3, ci_1, ci_2, ci_3, cf_1, cf_2, cf_3;
  int chi_i, chi_f, eps_f, eps_i, baryons;
  int chi_in[4][4][4], chi_out[4][4][4], eps[3][3][3];
  Real factor, iso_fac = 0;
  Real *prop = malloc(nt * sizeof(*prop));
  complex psi_sq, diquark, diquark_temp;
  static char *bar_kind[3] = {"PROTON", "DELTA", "DELTA_N"};

  // Set up epsilon tensor
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      for (k = 0; k < 3; k++)
        eps[i][j][k]= 0;
    }
  }
  eps[0][1][2]= 1;
  eps[1][2][0]= 1;
  eps[2][0][1]= 1;
  eps[0][2][1]= -1;
  eps[1][0][2]= -1;
  eps[2][1][0]= -1;

  for (baryons = PROTON; baryons <= DELTA_N; baryons++) {
    // Initialize chi factors
    for (i = 0; i < 4; i++) {
      for (j = 0; j < 4; j++) {
        for (k = 0; k < 4; k++) {
          chi_in[i][j][k] = 0;
          chi_out[i][j][k] = 0;
        }
      }
    }
    switch(baryons) {
      case PROTON:
        chi_in[1][0][ispinp1] =  1;
        chi_in[0][1][ispinp1] = -1;
        chi_in[3][2][ispinp1] =  1;
        chi_in[2][3][ispinp1] = -1;

        chi_out[1][0][ispinp2] =  1;
        chi_out[0][1][ispinp2] = -1;
        chi_out[3][2][ispinp2] =  1;
        chi_out[2][3][ispinp2] = -1;

        iso_fac = 1.0;
        break;

      case DELTA:
        chi_in[1][0][ispinp1] = 1;
        chi_in[0][1][ispinp1] = 1;
        chi_in[3][2][ispinp1] = 1;
        chi_in[2][3][ispinp1] = 1;

        chi_out[1][0][ispinp2] = 1;
        chi_out[0][1][ispinp2] = 1;
        chi_out[3][2][ispinp2] = 1;
        chi_out[2][3][ispinp2] = 1;

        iso_fac = 1.0;
        break;

      case DELTA_N:
        chi_in[0][0][ispinp1] = 1;
        chi_in[2][2][ispinp1] = 1;

        chi_out[0][0][ispinp2] = 1;
        chi_out[2][2][ispinp2] = 1;

        iso_fac = 2.0;
        break;
    }

    // Clear baryon propagator
    for (t = 0; t < nt; t++)
      *(prop + t) = 0.0;

    /* begin sum over source quark spin and color (labelled by 'f' suffix) */
    FORALLSITES(isite, s) {
      my_t = s->t;

      sf_3 =  ispinp2;
      for (cf_3 = 0; cf_3 < Nc; cf_3++) {
        // Sum over source diquark components
        for (sf_1 = 0; sf_1 < Ns; sf_1++)
          for (sf_2 = 0; sf_2 < Ns; sf_2++)
          {
            chi_f = chi_out[sf_1][sf_2][sf_3];
            if (chi_f != 0)
            {
              for (cf_1 = 0; cf_1 < Nc; cf_1++)
                for (cf_2 = 0; cf_2 < Nc; cf_2++)
                {
                  eps_f = eps[cf_1][cf_2][cf_3];
                  if (eps_f != 0)
                  {


                    /*  combine the sink colors and spins of these three propagators
                        to compose the outgoing baryon */

                    /* Sum over sink quark color and spin (labelled by 'i' suffix)*/
                    for (si_3 = 0; si_3 < Ns; si_3++)
                      for (ci_3 = 0; ci_3 < Nc; ci_3++)
                      {
                        /* begin sum over sink diquark components */
                        for (si_1 = 0; si_1 < Ns; si_1++)
                          for (si_2 = 0; si_2 < Ns; si_2++)
                          {
                            chi_i = chi_in[si_1][si_2][si_3];
                            if (chi_i != 0)
                            {
                              for (ci_1 = 0; ci_1 < Nc; ci_1++)
                                for (ci_2 = 0; ci_2 < Nc; ci_2++)
                                {
                                  eps_i = eps[ci_1][ci_2][ci_3];
                                  if (eps_i != 0) {
                                    factor = (Real)(eps_f * eps_i * chi_i * chi_f);

                                    /*  build the 2-3 diquark direct term */
                                    CMUL(
                                        ((wilson_matrix *)F_PT(s, src2))->d[sf_2].c[cf_2].d[si_2].c[ci_2],
                                        ((wilson_matrix *)F_PT(s, src3))->d[sf_3].c[cf_3].d[si_3].c[ci_3],
                                        diquark_temp);
                                    CMULREAL(diquark_temp, factor, diquark);

                                    /*  next, build the 2-3 diquark exchange term and subtract it from the
                                        diquark */
                                    CMUL(
                                        ((wilson_matrix *)F_PT(s, src2))->d[sf_2].c[cf_2].d[si_3].c[ci_3],
                                        ((wilson_matrix *)F_PT(s, src3))->d[sf_3].c[cf_3].d[si_2].c[ci_2],
                                        diquark_temp);
                                    factor = factor * iso_fac;
                                    CMULREAL(diquark_temp, factor, diquark_temp);
                                    CSUB(diquark, diquark_temp, diquark);

                                    /*  finally tie the diquark to quark 1 to form the baryon */
                                    CMUL(
                                        ((wilson_matrix *)F_PT(s, src1))->d[sf_1].c[cf_1].d[si_1].c[ci_1],
                                        diquark, psi_sq);
                                    /* and accumulate into the baryon propagator */
                                    *(prop  + my_t) += psi_sq.real;

                                  }/* end if eps_i */
                                } /* Close ci_1, ci_2 */
                            } /* end if chi_i */
                          } /* Close si_1, si_2 */
                      } /* Close  si_3, ci_3 */
                    /*  this ends all loops over the sink */

                  } /* end if eps_f */
                }  /* end do cf_1, cf_2 */
            }  /*end if chi_f */
          }  /*end do sf_1, sf_2 */
      }  /*end do cf_3   (and sf_3, if used) */
    }  /*end loop for sites */

    // Print the baryon propagators
    for (t = 0; t < nt; t++) {
      g_floatsum(prop + t);
      node0_printf("%s %d %.8g\n",
                   bar_kind[baryons], t, (double)(*(prop + t)));
    }
  }
  free(prop);
}
#undef Nc
#undef Ns
#undef PROTON
#undef DELTA
#undef DELTA_N
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int spectrum() {
  register int i;
  register site *s;
  int iters = 0, spin, color;

  for (spin = 0; spin < 4; spin++) {
    for (color = 0; color < NCOL; color++) {
#ifdef DEBUG_CHECK
      node0_printf("spectrum spin = %d, color = %d\n", spin, color);
#endif
      source(chi[1], color, spin, WALL, 0, 0, 0, 0, 0.25);

#ifdef DEBUG_CHECK
      printf("Dumping source...\n");
      FORALLSITES(i, s) {
        printf("chi[1](%d, %d, %d, %d):\n", s->x, s->y, s->z, s->t);
        dump_wvec(&(s->chi[1]));
      }

      printf("Dumping link_1...\n");
      FORALLSITES(i, s) {
        printf("link_1(%d, %d, %d, %d):\n", s->x, s->y, s->z, s->t);
        dumpmat(&(s->link[1]));
      }
#endif

      // chi[0] <-- Mdag.chi[1]
      make_clov(CKU0);
      make_clovinv(ODD);
      fermion_op(chi[1], chi[0], MINUS, EVEN);

      // Invert Mdag.M with zero initial guess and result in psi[0]
      FOREVENSITES(i, s)
        clear_wvec(&(psi[0][i]));
      iters += congrad(0, 0.0, EVEN);

      // Repeat the steps above for odd sites
      free_clov();
      make_clov(CKU0);
      make_clovinv(EVEN);
      fermion_op(chi[1], chi[0], MINUS, ODD);
      FORODDSITES(i, s)
        clear_wvec(&(psi[0][i]));
      iters += congrad(0, 0.0, ODD);
      free_clov();

#ifdef DEBUG_CHECK
      printf("Dumping Mdag.source...\n");
      FORALLSITES(i, s) {
        printf("chi[0](%d, %d, %d, %d):\n", s->x, s->y, s->z, s->t);
        dump_wvec(&(chi[0][i]));
      }

      printf("Dumping propagator...\n");
      FORALLSITES(i, s) {
        printf("psi[0](%d, %d, %d, %d):\n", s->x, s->y, s->z, s->t);
        dump_wvec(&(psi[0][i]));
      }
#endif

      // Copy result to quark_propagator
      FORALLSITES(i, s)
        copy_wvec(&(psi[0][i]), &(s->quark_propagator.d[spin].c[color]));

      // TODO: The 'DSL' mesons change far more than the other results
      //       when using congrad rather than wilson_invert...
      //       Comment out for now
//#if 0
      /* For extra contributions needed for clover "rotated" fields */
      /* contributing to the pion propagator--gamma5 gamma5 correlator and
         gamma5 gamma5.gamma0 correlator */
      /* These refinements are required for the axial Ward identity */

      /* These extra contributions are of the form
         psibar <-dslash Gammas psi  and psibar Gammas dslash->  psi  */
      /* Where dslash is the naive operator (without the Wilson term) */

      // dslash on the propagator
      dslash(psi[0], mp, PLUS, EVENANDODD);
      dslash(psi[0], tempwvec, MINUS, EVENANDODD);
      // From subtraction we get 2Dslash
      FORALLSITES(i, s) {
        sub_wvec(&(mp[i]), &(tempwvec[i]),
                 &s->rotated_propagator.d[spin].c[color]);
      }
//#endif
    }
  }
  meson(F_OFFSET(quark_propagator), F_OFFSET(quark_propagator), "ONE_");
//  meson(F_OFFSET(quark_propagator), F_OFFSET(rotated_propagator), "DSL_");
  baryon(F_OFFSET(quark_propagator), F_OFFSET(quark_propagator),
         F_OFFSET(quark_propagator));

  return iters;
}
// -----------------------------------------------------------------
