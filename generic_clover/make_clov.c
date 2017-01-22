// -----------------------------------------------------------------
// Prepare the clover term
#include "generic_clover_includes.h"
#include "../include/prefetch.h"
#define FETCH_UP 1
#define LOOPEND
#include "../include/loopend.h"

// Use global clover term when working with only one clover term at a time
static clover *global_clov;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Allocate space for a clover term
clover *create_clov(void) {
  clover *my_clov = malloc(sizeof(*my_clov));
  my_clov->clov = malloc(sites_on_node * sizeof(triangular));
  my_clov->clov_diag = malloc(sites_on_node * sizeof(diagonal));

  if (my_clov->clov == NULL || my_clov->clov_diag == NULL) {
    printf("create_clov(%d): malloc failed\n", this_node);
    free(my_clov);
    return NULL;
  }
  else
    return my_clov;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute the clover field
void compute_clov(clover *my_clov, Real Clov_c) {
  register int i, j, k, jk, jk2;
  register site *s;
  register complex tc;
  triangular *clov = my_clov->clov;
  diagonal *clov_diag = my_clov->clov_diag;
  su3_matrix *f_mn = malloc(sites_on_node * sizeof(*f_mn));
  char myname[] = "make_clov";

  if (f_mn == NULL) {
    printf("%s(%d): Can't malloc f_mn\n", myname, this_node);
    terminate(1);
  }

/* The clover term
*                                           (X | 0)
* sum_{mu < nu} sigma_{mu, nu} F_{mu, nu} = (-----)
*                                           (0 | Y)
*
* is block diagonal in the Weyl basis used, and hermitian. Together with the 1|
* we will store it in a lower complex triangular (block) matrix (without
* diagonal) and a real diagonal. The blocks go over Dirac indices 0, 1 and 2, 3.
*
*   With
*
*   sigma_{mu, nu} = i/2 [gamma_{mu}, gamma_{nu}]
*
*   and with our gamma matrices (see libraries/mb_gamma.c)
*
*   1234 = XYZT    and   Pauli sigma_1, sigma_2, sigma_3
*
*   sigma_23 =  diag(sigma_1,  sigma_1)
*   sigma_31 =  diag(-sigma_2, -sigma_2)
*   sigma_12 =  diag(sigma_3,  sigma_3)
*   sigma_14 =  diag(-sigma_1,  sigma_1)
*   sigma_24 =  diag(sigma_2, -sigma_2)
*   sigma_34 =  diag(-sigma_3,  sigma_3)
*
*
* Here X = sigma_1 (F_23 - F_14) + sigma_2 (F_13 + F_24) + sigma_3 (F_12 - F_34)
* and  Y = sigma_1 (F_23 + F_14) + sigma_2 (F_13 - F_24) + sigma_3 (F_12 + F_34)
*
* This has to be multiplied by the "clover coefficient" and subtracted from 1|.
*
* Note that F_mn = f_mn/i = -i f_mn!
*
* (The indices above start from 1, in the program they start from 0,
* i.e. subtract -1 from the above) */

/* The clover term R is a hermitian 12 x 12 matrix in color and Dirac spin.
   If we use ab for color and ij for spin, with a = (0, 1, 2) and
   i = (0, 1, 2, 3), then we store the matrix R(ai, bj) as follows with
   the color index varying most rapidly:

   Upper left block, showing blocks of SU(3) matrices

   [ 1 + ci(f_01 - f_23)                ci(f_12 - f_03) + c(f_02 + f_13) ]
   [ ci(f_12 - f_03) - c(f_02 + f_13)   1 - ci(f_01 - f_23)              ]

   Lower right block, showing blocks of SU(3) matrices

   [ 1 + ci(f_01 + f_23)                 ci(f_12 + f_03) + c(f_02 - f_13)  ]
   [ ci(f_12 + f_03) - c(f_02 - f_13)    1 - ci(f_01 + f_23)               ]

   This hermitian matrix is decomposed into a lower triangular matrix
   and its diagonals.  The diagonals are real.

   clov_diag[site].di[0][0-5] holds the diagonals from the upper left
   clov_diag[site].di[1][0-5] holds the diagonals from the lower right

   clov_diag[site].di[0][a]     = R(a0, a0) = 1 - c(Im f_01 - Im f_23)(a, a)
   clov_diag[site].di[0][a + 3] = R(a1, a1) = 1 + c(Im f_01 - Im f_23)(a, a)
   clov_diag[site].di[1][a]     = R(a2, a2) = 1 - c(Im f_01 + Im f_23)(a, a)
   clov_diag[site].di[1][a + 3] = R(a3, a3) = 1 + c(Im f_01 + Im f_23)(a, a)

   clov[site].tr[0][jk]:
      Holds the lower triangular elements of the upper left block in
      the pattern

        jk = ...
   [          |        ]      [                      |                     ]
   [  0       |        ]      [  1 + ci(f_01 - f_23) | ci(f_12 - f_03)     ]
   [  1  2    |        ]      [                      |   + c(f_02 + f_13)  ]
   ---------------------   =  ----------------------------------------------
   [  3  4  5 |        ]      [                      |                     ]
   [  6  7  8 |  9     ]      [ ci(f_12 - f_03)      | 1 - ci(f_01 - f_23) ]
   [ 10 11 12 | 13 14  ]      [    - c(f_02 + f_13)  |                     ]

   clov[site].tr[1][jk]:
      Holds the lower triangular elements of the lower right block in
      the pattern

        jk = ...
   [          |        ]      [                      |                     ]
   [  0       |        ]      [  1 + ci(f_01 + f_23) | ci(f_12 + f_03)     ]
   [  1  2    |        ]      [                      |   + c(f_02 - f_13)  ]
   ---------------------   =  ----------------------------------------------
   [  3  4  5 |        ]      [                      |                     ]
   [  6  7  8 |  9     ]      [ ci(f_12 + f_03)      | 1 - ci(f_01 + f_23) ]
   [ 10 11 12 | 13 14  ]      [    - c(f_02 - f_13)  |                     ]
   */

  f_mu_nu(f_mn, 0, 1);
  FORALLSITES(i, s)
    scalar_mult_su3_matrix(&(f_mn[i]), Clov_c, &(f_mn[i]));

  jk = 0;
  for (j = 0; j < DIMF; j++) {
    FORALLSITES(i, s) {
      clov_diag[i].di[0][j] = 1.0 - f_mn[i].e[j][j].imag;
      clov_diag[i].di[0][j + DIMF] = 1.0 + f_mn[i].e[j][j].imag;
      clov_diag[i].di[1][j] = 1.0 - f_mn[i].e[j][j].imag;
      clov_diag[i].di[1][j + DIMF] = 1.0 + f_mn[i].e[j][j].imag;
    }
    jk2 = (j + DIMF) * (j + DIMF-1)/2 + DIMF;
    for (k = 0; k < j; k++) {
      FORALLSITES(i, s) {
        CMUL_I(f_mn[i].e[j][k], clov[i].tr[0][jk]);
        CMUL_MINUS_I(f_mn[i].e[j][k], clov[i].tr[0][jk2]);
        CMUL_I(f_mn[i].e[j][k], clov[i].tr[1][jk]);
        CMUL_MINUS_I(f_mn[i].e[j][k], clov[i].tr[1][jk2]);
      }
      jk++;
      jk2++;
    }
  }

  f_mu_nu(f_mn, 2, 3);
  FORALLSITES(i, s)
    scalar_mult_su3_matrix(&(f_mn[i]), Clov_c, &(f_mn[i]));

  jk = 0;
  for (j = 0; j < DIMF; j++) {
    FORALLSITES(i, s) {
      clov_diag[i].di[0][j] += f_mn[i].e[j][j].imag;
      clov_diag[i].di[0][j + DIMF] -= f_mn[i].e[j][j].imag;
      clov_diag[i].di[1][j] -= f_mn[i].e[j][j].imag;
      clov_diag[i].di[1][j + DIMF] += f_mn[i].e[j][j].imag;
    }
    jk2 = (j + DIMF)*(j + DIMF-1)/2 + DIMF;
    for (k = 0; k < j; k++) {
      FORALLSITES(i, s) {
        CMUL_MINUS_I(f_mn[i].e[j][k], tc);
        CSUM(clov[i].tr[0][jk], tc);
        CDIF(clov[i].tr[0][jk2], tc);
        CDIF(clov[i].tr[1][jk], tc);
        CSUM(clov[i].tr[1][jk2], tc);
      }
      jk++;
      jk2++;
    }
  }

  f_mu_nu(f_mn, 1, 2);
  FORALLSITES(i, s)
    scalar_mult_su3_matrix(&(f_mn[i]), Clov_c, &(f_mn[i]));

  for (j = 0; j < DIMF; j++) {
    jk = (j + DIMF)*(j + DIMF-1)/2;
    for (k = 0; k < DIMF; k++) {
      FORALLSITES(i, s) {
        CMUL_I(f_mn[i].e[j][k], clov[i].tr[0][jk]);
        CMUL_I(f_mn[i].e[j][k], clov[i].tr[1][jk]);
      }
      jk++;
    }
  }

  f_mu_nu(f_mn, 0, 3);
  FORALLSITES(i, s)
    scalar_mult_su3_matrix(&(f_mn[i]), Clov_c, &(f_mn[i]));

  for (j = 0; j < DIMF; j++) {
    jk = (j + DIMF)*(j + DIMF-1)/2;
    for (k = 0; k < DIMF; k++) {
      FORALLSITES(i, s) {
        CMUL_MINUS_I(f_mn[i].e[j][k], tc);
        CSUM(clov[i].tr[0][jk], tc);
        CDIF(clov[i].tr[1][jk], tc);
      }
      jk++;
    }
  }

  f_mu_nu(f_mn, 0, 2);
  FORALLSITES(i, s)
    scalar_mult_su3_matrix(&(f_mn[i]), Clov_c, &(f_mn[i]));

  for (j = 0; j < DIMF; j++) {
    jk = (j + DIMF) * (j + DIMF - 1) / 2;
    for (k = 0; k < DIMF; k++) {
      FORALLSITES(i, s) {
        CDIF(clov[i].tr[0][jk], f_mn[i].e[j][k]);
        CDIF(clov[i].tr[1][jk], f_mn[i].e[j][k]);
      }
      jk++;
    }
  }

  f_mu_nu(f_mn, 1, 3);
  FORALLSITES(i, s)
    scalar_mult_su3_matrix(&(f_mn[i]), Clov_c, &(f_mn[i]));

  for (j = 0; j < DIMF; j++) {
    jk = (j + DIMF)*(j + DIMF-1)/2;
    for (k = 0; k < DIMF; k++) {
      FORALLSITES(i, s) {
        CDIF(clov[i].tr[0][jk], f_mn[i].e[j][k]);
        CSUM(clov[i].tr[1][jk], f_mn[i].e[j][k]);
      }
      jk++;
    }
  }
  free(f_mn);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Invert the clover term on a specified sublattice
// Return tr(log(A)) for HMC
double compute_clovinv(clover *my_clov, int parity) {
  triangular *clov = my_clov->clov;
  diagonal *clov_diag = my_clov->clov_diag;
  register int i, b, j, k, kj, l, lk, lj, jl;
  register site *s;
  Real f_diag[2 * DIMF];
  double trlogA = 0.0;
  complex v1[2 * DIMF], tc, sum;

  // Take the inverse on the specified sublattice for each of the 2 blocks
  FORSOMEPARITY(i, s, parity) {
    for (b = 0; b < 2; b++) {
      // Cholesky decompose
      // Do in place with temporary diagonal.
      // Algorithm 4.2.2 of "Matrix Computations" (Golub & Van Loan, 1989)
      for (j = 0; j < 2 * DIMF; j++) {

        clov_diag[i].di[b][j] = sqrt((double)(clov_diag[i].di[b][j]));
        f_diag[j] = 1.0 / (clov_diag[i].di[b][j]);
        for (k = j + 1; k < 2 * DIMF; k++) {
          kj = k * (k - 1) / 2 + j;
          CMULREAL(clov[i].tr[b][kj], f_diag[j], clov[i].tr[b][kj]);
        }

        for (k = j + 1; k < 2 * DIMF; k++) {
          kj = k*(k-1)/2 + j;
          CMUL_J(clov[i].tr[b][kj], clov[i].tr[b][kj], tc);
          clov_diag[i].di[b][k] -= tc.real;
          for (l = k + 1; l < 2 * DIMF; l++) {
            lj = l*(l-1)/2 + j;
            lk = l*(l-1)/2 + k;
            CMUL_JDIF(clov[i].tr[b][lj], clov[i].tr[b][kj], clov[i].tr[b][lk]);
          }
        }
      }

      /* Accumulate trlogA */
      for (j = 0; j < 2 * DIMF; j++)
        trlogA += (double)2.0 * log((double)(clov_diag[i].di[b][j]));

      /* Now use forward and backward substitution to construct inverse */
      for (k = 0; k < 2 * DIMF; k++) {
        for (l = 0; l < k; l++)
          v1[l] = cmplx(0.0, 0.0);

        /* Forward substitute */
        v1[k] = cmplx(f_diag[k], 0.0);
        for (l = k + 1; l < 2 * DIMF; l++) {
          sum = cmplx(0.0, 0.0);
          for (j = k; j < l; j++) {
            lj = l * (l - 1) / 2 + j;
            CMULDIF(clov[i].tr[b][lj], v1[j], sum);
          }
          CMULREAL(sum, f_diag[l], v1[l]);
        }

        /* Backward substitute */
        l = 2 * DIMF - 1;
        CMULREAL(v1[l], f_diag[l], v1[l]);
        for (l = 2 * DIMF - 2; l >= k; l--) {
          sum = v1[l];
          for (j = l + 1; j < 2 * DIMF; j++) {
            jl = j * (j - 1) / 2 + l;
            CMULJ_DIF(clov[i].tr[b][jl], v1[j], sum);
          }
          CMULREAL(sum, f_diag[l], v1[l]);
        }

        /* Overwrite column k */
        clov_diag[i].di[b][k] = v1[k].real;
        for (l = k + 1; l < 2 * DIMF; l++) {
          lk = l * (l - 1) / 2 + k;
          clov[i].tr[b][lk] = v1[l];
        }
      }
    }
  } END_LOOP

  g_doublesum(&trlogA);
  return trlogA;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
/* Multiply a Wilson vector (spin, color) with a block-diagonal hermition
   matrix stored as a complex lower triangular matrix (without diagonal)
   and a real diagonal. The blocks are spin indices 0, 1 and 2, 3. */

/* To simplify the code with the above block structure, introduce
   the somewhat dirty structure, equivalent to a wilson_vector: */
typedef struct { complex b[2][2 * DIMF]; } wilson_block_vector;

// Recast src and dest as wilson_block_vector
void mult_this_ldu(clover *my_clov, wilson_vector *src,
                   wilson_vector *dest, int parity) {

  register int i, b, j, k, jk, kj;
  register site *s;
  triangular *clov = my_clov->clov;
  diagonal *clov_diag = my_clov->clov_diag;
  wilson_block_vector *srcb = (wilson_block_vector *)src;
  wilson_block_vector *destb = (wilson_block_vector *)dest;

  FORSOMEPARITY(i, s, parity) {
    if (i < loopend - FETCH_UP) {
      prefetch_W(&srcb[i + FETCH_UP]);
      prefetch_W(&destb[i + FETCH_UP]);
    }
    for (b = 0; b < 2; b++) {
      for (j = 0; j < 2*DIMF; j++) {
        /* diagonal part */
        CMULREAL(srcb[i].b[b][j], clov_diag[i].di[b][j], destb[i].b[b][j]);

        /* lower triangular part */
        jk = j * (j - 1) / 2;
        for (k = 0; k < j; k++) {
          CMULSUM(clov[i].tr[b][jk], srcb[i].b[b][k], destb[i].b[b][j]);
          jk++;
        }

        /* upper triangular part */
        for (k = j + 1; k < 2 * DIMF; k++) {
          kj = k * (k - 1) / 2 + j;
          CMULJ_SUM(clov[i].tr[b][kj], srcb[i].b[b][k], destb[i].b[b][j]);
        }
      }
    }
  } END_LOOP
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Left-multiply a color-Dirac matrix stored in ldu form
// by sigma_mu_nu = gamma_mu gamma_nu, where
/*
sigma(XUP, YUP)          sigma(XUP, ZUP)          sigma(XUP, TUP)
        -i  0  0  0              0 -1  0  0              0  i  0  0
         0  i  0  0              1  0  0  0              i  0  0  0
         0  0 -i  0              0  0  0 -1              0  0  0 -i
         0  0  0  i              0  0  1  0              0  0 -i  0

sigma(YUP, ZUP)          sigma(YUP, TUP)          sigma(ZUP, TUP)
         0 -i  0  0              0 -1  0  0              i  0  0  0
        -i  0  0  0              1  0  0  0              0 -i  0  0
         0  0  0 -i              0  0  0  1              0  0 -i  0
         0  0 -i  0              0  0 -1  0              0  0  0  i

and sigma(nu, mu) = -sigma(mu, nu), and sums over the dirac indices.
*/
// Used in clover_diagonal

/* triang, diag are input & contain the color-dirac matrix
   mat is output: the resulting su3_matrix  */
// Put result in tempmat
void tr_sigma_ldu_mu_nu(int mu, int nu) {
  register int i, j, k, jk, jk2, kj;
  register int mm = 0, nn = 0;  /* dummy directions */
  register Real pm = 0;
  register site *s;
  Real tr;
  complex tc, tc2;
  triangular *clov = global_clov->clov;
  diagonal *clov_diag = global_clov->clov_diag;
  char myname[] = "tr_sigma_ldu_mu_nu";

  /* take care of the case mu > nu by flipping them and mult. by -1 */
  if (mu < nu) {
    mm = mu;
    nn = nu;
    pm = 1.0;
  }
  else if (mu > nu) {
    mm = nu;
    nn = mu;
    pm = -1.0;
  }
  else
    printf("BAD CALL by %s: mu =%d, nu =%d\n", myname, mu, nu);

  switch(mm) {
    case(XUP):
      switch(nn) {
        case(YUP):
          FORODDSITES(i, s) {
            /* diagonal part */
            for (j = 0; j < DIMF; j++) {
              tr = clov_diag[i].di[0][j + DIMF] - clov_diag[i].di[0][j];
              tr -= clov_diag[i].di[1][j];
              tr += clov_diag[i].di[1][j + DIMF];
              tempmat[i].e[j][j].real = 0.0;
              tempmat[i].e[j][j].imag = tr;
            }
            /* triangular part */
            jk = 0;
            for (j = 1; j < DIMF; j++) {
              jk2 = (j + DIMF) * (j + DIMF - 1) / 2 + DIMF;
              for (k = 0; k < j; k++) {
                CSUB(clov[i].tr[0][jk2], clov[i].tr[0][jk], tc);
                CDIF(tc, clov[i].tr[1][jk]);
                CSUM(tc, clov[i].tr[1][jk2]);
                CMUL_I(tc, tempmat[i].e[j][k]);
                CONJG(tc, tc);
                CMUL_I(tc, tempmat[i].e[k][j]);
                jk++;
                jk2++;
              }
            }
          }
          break;
        case(ZUP):
          FORODDSITES(i, s) {
            for (j = 0; j < DIMF; j++) {
              jk = (j + DIMF)*(j + DIMF-1)/2;
              for (k = 0; k < DIMF; k++) {
                kj = (k + DIMF)*(k + DIMF-1)/2 + j;
                CONJG(clov[i].tr[0][kj], tc);
                CDIF(tc, clov[i].tr[0][jk]);
                CONJG(clov[i].tr[1][kj], tc2);
                CSUM(tc, tc2);
                CSUB(tc, clov[i].tr[1][jk], tempmat[i].e[j][k]);
                jk++;
              }
            }
          }
          break;
        case(TUP):
          FORODDSITES(i, s) {
            for (j = 0; j < DIMF; j++) {
              jk = (j + DIMF)*(j + DIMF-1)/2;
              for (k = 0; k < DIMF; k++) {
                kj = (k + DIMF)*(k + DIMF-1)/2 + j;
                CONJG(clov[i].tr[0][kj], tc);
                CSUM(tc, clov[i].tr[0][jk]);
                CONJG(clov[i].tr[1][kj], tc2);
                CDIF(tc, tc2);
                CDIF(tc, clov[i].tr[1][jk]);
                CMUL_I(tc, tempmat[i].e[j][k]);
                jk++;
              }
            }
          }
          break;
        default:
          printf("BAD CALL in %s: mu =%d, nu =%d\n", myname, mu, nu);
      }
      break;
    case(YUP):
      switch(nn) {
        case(ZUP):
          FORODDSITES(i, s) {
            for (j = 0; j < DIMF; j++) {
              jk = (j + DIMF)*(j + DIMF-1)/2;
              for (k = 0; k < DIMF; k++) {
                kj = (k + DIMF)*(k + DIMF-1)/2 + j;
                CONJG(clov[i].tr[0][kj], tc);
                CSUM(tc, clov[i].tr[0][jk]);
                CONJG(clov[i].tr[1][kj], tc2);
                CSUM(tc, tc2);
                CSUM(tc, clov[i].tr[1][jk]);
                CMUL_MINUS_I(tc, tempmat[i].e[j][k]);
                jk++;
              }
            }
          }
          break;
        case(TUP):
          FORODDSITES(i, s) {
            for (j = 0; j < DIMF; j++) {
              jk = (j + DIMF)*(j + DIMF-1)/2;
              for (k = 0; k < DIMF; k++) {
                kj = (k + DIMF)*(k + DIMF-1)/2 + j;
                CONJG(clov[i].tr[0][kj], tc);
                CDIF(tc, clov[i].tr[0][jk]);
                CONJG(clov[i].tr[1][kj], tc2);
                CDIF(tc, tc2);
                CADD(tc, clov[i].tr[1][jk], tempmat[i].e[j][k]);
                jk++;
              }
            }
          }
          break;
        default:
          printf("BAD CALL in %s: mu =%d, nu =%d\n", myname, mu, nu);
      }
      break;
    case(ZUP):
      switch(nn) {
        case(TUP):
          FORODDSITES(i, s) {
            /* diagonal part */
            for (j = 0; j < DIMF; j++) {
              tr = clov_diag[i].di[0][j] - clov_diag[i].di[0][j + DIMF];
              tr -= clov_diag[i].di[1][j];
              tr += clov_diag[i].di[1][j + DIMF];
              tempmat[i].e[j][j].real = 0.0;
              tempmat[i].e[j][j].imag = tr;
            }
            /* triangular part */
            jk = 0;
            for (j = 1; j < DIMF; j++) {
              jk2 = (j + DIMF)*(j + DIMF-1)/2 + DIMF;
              for (k = 0; k < j; k++) {
                CSUB(clov[i].tr[0][jk], clov[i].tr[0][jk2], tc);
                CDIF(tc, clov[i].tr[1][jk]);
                CSUM(tc, clov[i].tr[1][jk2]);
                CMUL_I(tc, tempmat[i].e[j][k]);
                CONJG(tc, tc);
                CMUL_I(tc, tempmat[i].e[k][j]);
                jk++;
                jk2++;
              }
            }
          }
          break;
        default:
          printf("BAD CALL in %s: mu =%d, nu =%d\n", myname, mu, nu);
      }
      break;
    default:
      printf("BAD CALL in %s: mu =%d, nu =%d\n", myname, mu, nu);
  }

  // Multiply by pm = +/- 1
  if (pm < 0) {
    FORODDSITES(i, s) {
      for (j = 0; j < DIMF; j++) {
        for (k = 0; k < DIMF; k++)
          CNEGATE(tempmat[i].e[j][k], tempmat[i].e[j][k]);
      }
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Make only the global clover term
// Compute R = 1 - i * CKU0 sigma_{\mu\nu} F_{\mu\nu} on each site
void make_clov(Real Clov_c) {
  global_clov = create_clov();
  if (global_clov == NULL)
    terminate(1);

  compute_clov(global_clov, Clov_c);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Invert only the global clover term
// "clov" and "clov_diag" now contain R^(-1) on parity sites,
// still with R itself on other_parity sites
double make_clovinv(int parity) {
  return compute_clovinv(global_clov, parity);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Multiply only by the global clover term
// src and dest will be recast as wilson_block_vector
void mult_ldu(wilson_vector *src, wilson_vector *dest, int parity) {
  mult_this_ldu(global_clov, src, dest, parity);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void free_this_clov(clover *my_clov) {
  free(my_clov->clov);
  free(my_clov->clov_diag);
  free(my_clov);
}

void free_clov() {
  free_this_clov(global_clov);
}
// -----------------------------------------------------------------
