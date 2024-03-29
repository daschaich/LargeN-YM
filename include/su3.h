// -----------------------------------------------------------------
// Defines and subroutine declarations for SU(N) gauge theory
// TODO: A bit of inline doc cleanup in progress...
#ifndef _SUN_H
#define _SUN_H

#include "../include/complex.h"
#include "../include/random.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Set number of colors
#define NCOL 3
#define N_OFFDIAG (NCOL * (NCOL - 1) / 2)

#ifdef FAST
  #if (NCOL != 3)
    #error "FAST only works if NCOL = 3"
  #endif
#endif

typedef struct { fcomplex e[NCOL][NCOL]; } fmatrix;
typedef struct { fcomplex c[NCOL]; } fvector;

// Anti-hermitian matrices for general NCOL
typedef struct {
  fcomplex m[N_OFFDIAG];
  float im_diag[NCOL];
} fanti_hermitmat;

typedef struct { dcomplex e[NCOL][NCOL]; } dmatrix;
typedef struct { dcomplex c[NCOL]; } dvector;
typedef struct {
  dcomplex m[N_OFFDIAG];
  double im_diag[NCOL];
} danti_hermitmat;

#if PRECISION == 1
#define matrix    fmatrix
#define vector    fvector
#define anti_hermitmat  fanti_hermitmat
#else
#define matrix    dmatrix
#define vector    dvector
#define anti_hermitmat  danti_hermitmat
#endif

// Need SU(2) matrices for any SU(N), e.g. for gauge hits when gauge-fixing
typedef struct { complex e[2][2]; } su2_matrix;

/*
* TODO: CLEANUP IN PROGRESS...
*
* ROUTINES FOR MATRIX OPERATIONS
* void make_anti_hermitian(matrix *m3,  anti_hermitmat *ah3)
* file "make_ahmat.c"
* void random_anti_hermitian(anti_hermitmat *mat_antihermit, double_prn *prn_pt)
* (prn_pt passed through to myrand())
* file "rand_ahmat.c"
* void uncompress_anti_hermitian(anti_hermitmat *mat_anti, matrix *mat)
* file "uncmp_ahmat.c"
* void compress_anti_hermitian(matrix *mat, anti_hermitmat *mat_anti)
* file "cmp_ahmat.c"
*
* complex det(matrix *a)
* file "det.c"
*
* ROUTINES MIXING SU(2) and SU(N)
* void left_su2_hit_n(su2_matrix *u, int p, int q, matrix *link)
*       file "l_su2_hit_n.c"
* void right_su2_hit_a(su2_matrix *u, int p, int q, matrix *link)
*       file "r_su2_hit_a.c"
* void dump_su2(su2_matrix *u)
*       file "dump_su2.c"
*/
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Subroutine definitions
// Vector operations
// In file clearvec.c
void clearvec(vector *v);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Anti-hermitian matrix routines
void make_anti_hermitian(matrix *m, anti_hermitmat *ah);
void random_anti_hermitian(anti_hermitmat *ah, double_prn *prn_pt);
void uncompress_anti_hermitian(anti_hermitmat *ah, matrix *m);
void compress_anti_hermitian(matrix *m, anti_hermitmat *ah);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Fundamental matrix operations
// In file dump_mat.c
void dump_mat(matrix *m);

// In file clear_mat.c
void clear_mat(matrix *m);

// In file trace.c
complex trace(matrix *a);
void trace_sum(matrix *a, complex *c);

// In file adjoint.c
void adjoint(matrix *a, matrix *b);

// In file s_a_d_mat.c
void scalar_add_diag(matrix *a, Real s);

// In file cs_a_d_mat.c
void c_scalar_add_diag(matrix *a, complex *s);

// In file tr_prod.c
Real realtrace(matrix *a, matrix *b);
Real realtrace_nn(matrix *a, matrix *b);
void realtrace_sum(matrix *a, matrix *b, Real *c);
complex complextrace(matrix *a, matrix *b);

// In file mat_copy.c
void mat_copy(matrix *a, matrix *b);

// In file addmat.c
void add_mat(matrix *a, matrix *b, matrix *c);
void sum_mat(matrix *b, matrix *c);
void sub_mat(matrix *a, matrix *b, matrix *c);
void dif_mat(matrix *b, matrix *c);

// In file s_m_mat.c
void scalar_mult_mat(matrix *src, Real scalar, matrix *dest);
void scalar_mult_add_mat(matrix *a, matrix *b, Real s, matrix *c);
void scalar_mult_sum_mat(matrix *b, Real s, matrix *c);
void scalar_mult_dif_mat(matrix *b, Real s, matrix *c);

// In file cs_m_mat.c
void c_scalar_mult_mat(matrix *b, complex *s, matrix *c);
void c_scalar_mult_sum_mat(matrix *b, complex *s, matrix *c);
void c_scalar_mult_add_mat(matrix *a, matrix *b, complex *s,
                           matrix *c);

// In file m_mat_nn.c
void mult_nn(matrix *a, matrix *b, matrix *c);
void mult_nn_sum(matrix *a, matrix *b, matrix *c);

// In file m_mat_na.c
void mult_na_sum(matrix *a, matrix *b, matrix *c);
void mult_na(matrix *a, matrix *b, matrix *c);

// In file m_mat_an.c
void mult_an(matrix *a, matrix *b, matrix *c);
void mult_an_sum(matrix *a, matrix *b, matrix *c);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Routines mixing SU(2) and U(N)
// In file m_su2_mat_vec_n.c
void mult_su2_mat_vec_elem_n(su2_matrix *u, complex *x0, complex *x1);

// In file m_su2_mat_vec_a.c
void mult_su2_mat_vec_elem_a(su2_matrix *u, complex *x0, complex *x1);

void left_su2_hit_n(su2_matrix *u, int p, int q, matrix *link);
void right_su2_hit_a(su2_matrix *u, int p, int q, matrix *link);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Miscellaneous routines
complex det(matrix *a);

void dump_su2(su2_matrix *u);

// In file gaussrand.c
Real gaussian_rand_no(double_prn *prn_pt);

// In file byterevn.c
#include "../include/int32type.h"
void byterevn(int32type w[], int n);
void byterevn64(int32type w[], int n);
#endif
// -----------------------------------------------------------------
