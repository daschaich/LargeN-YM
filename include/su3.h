// -----------------------------------------------------------------
// Defines and subroutine declarations for SU(N) gauge theory
// with various dimension-DIMF fermion reps
// The original names now refer to objects of dimension DIMF or DIMFxDIMF
// New objects with suffix _f have dimension NCOL or NCOLxNCOL
#ifndef _SUN_H
#define _SUN_H

#include "../include/complex.h"
#include "../include/random.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Irrep stuff
// I'm having some trouble #defining FREP as <text>
// Let's make an integer code for each irrep, as in macros.h and io_lat.h
#define FUNDAMENTAL    90
#define SYMMETRIC2     91
#define ANTISYMMETRIC2 92
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// FREP can take the values: FUNDAMENTAL, SYMMETRIC2, ANTISYMMETRIC2
#define NCOL 3
#define DIMF 3
#define FREP FUNDAMENTAL
#define N_OFFDIAG (NCOL * (NCOL - 1) / 2)

// Only have anti_hermitmat and explicit loops in libraries for N <= 4
#if NCOL > 4
  #error "NCOL > 4 not yet implemented!"
#endif

#if (NCOL != 3 || DIMF != 3)
  #ifdef FAST
    #error "FAST only works if NCOL = DIMF = 3"
  #endif
#endif

typedef struct { fcomplex e[NCOL][NCOL]; } fmatrix_f;
typedef struct { fcomplex e[DIMF][DIMF]; } fmatrix;
typedef struct { fcomplex c[NCOL]; } fvector_f;
typedef struct { fcomplex c[DIMF]; } fvector;

// Anti-hermitian matrices for general NCOL
typedef struct {
  fcomplex m[N_OFFDIAG];
  float im_diag[NCOL];
} fanti_hermitmat;

typedef struct { dcomplex e[NCOL][NCOL]; } dmatrix_f;
typedef struct { dcomplex e[DIMF][DIMF]; } dmatrix;
typedef struct { dcomplex c[NCOL]; } dvector_f;
typedef struct { dcomplex c[DIMF]; } dvector;
typedef struct {
  dcomplex m[N_OFFDIAG];
  double im_diag[NCOL];
} danti_hermitmat;

#if PRECISION == 1
#define matrix_f   fmatrix_f
#define matrix     fmatrix
#define vector_f   fvector_f
#define vector     fvector
#define anti_hermitmat fanti_hermitmat
#else
#define matrix_f   dmatrix_f
#define matrix     dmatrix
#define vector_f   dvector_f
#define vector     dvector
#define anti_hermitmat danti_hermitmat
#endif

/* SU(2) */
typedef struct { complex e[2][2]; } su2_matrix;

/* Wilson vectors */
/* e.g.                */
/* wilson_propagator prop;                           */
/* prop.c[ci].d[si].d[sf].c[cf]                      */
/* ----------------------->    complex               */
/* ----------------->          vector            */
/* ----------->                wilson_vector         */
/* ----->                      spin_wilson_vector    */
/* e.g.                */
/* wilson_matrix matr;                               */
/* matr.d[si].c[ci].d[sf].c[cf]                      */
/* ----------------------->    complex               */
/* ----------------->          vector            */
/* ----------->                wilson_vector         */
/* ----->                      color_wilson_vector   */

/* Object with two Dirac and two color indices. A given element
   of a "wilson_propagator" is accessed by
   object.c[color1].d[spin1].d[spin2].c[color2].real , etc.
   As alway, "d" denotes a Dirac index and "c" a color index.
   "1" refers to the source, "2" to the sink.

   Color indices here always run to DIMF, not NCOL
*/

typedef struct { fvector d[4]; } fwilson_vector;
typedef struct { fvector h[2]; } fhalf_wilson_vector;
typedef struct { fwilson_vector c[DIMF]; } fcolor_wilson_vector;
typedef struct { fwilson_vector d[4]; } fspin_wilson_vector;
typedef struct { fcolor_wilson_vector d[4]; } fwilson_matrix;
typedef struct { fspin_wilson_vector c[DIMF]; } fwilson_propagator;

typedef struct { dvector d[4]; } dwilson_vector;
typedef struct { dvector h[2]; } dhalf_wilson_vector;
typedef struct { dwilson_vector c[DIMF]; } dcolor_wilson_vector;
typedef struct { dwilson_vector d[4]; } dspin_wilson_vector;
typedef struct { dcolor_wilson_vector d[4]; } dwilson_matrix;
typedef struct { dspin_wilson_vector c[DIMF]; } dwilson_propagator;

#if PRECISION == 1
#define wilson_vector       fwilson_vector
#define half_wilson_vector  fhalf_wilson_vector
#define color_wilson_vector fcolor_wilson_vector
#define spin_wilson_vector  fspin_wilson_vector
#define wilson_matrix       fwilson_matrix
#define wilson_propagator   fwilson_propagator
#else
#define wilson_vector       dwilson_vector
#define half_wilson_vector  dhalf_wilson_vector
#define color_wilson_vector dcolor_wilson_vector
#define spin_wilson_vector  dspin_wilson_vector
#define wilson_matrix       dwilson_matrix
#define wilson_propagator   dwilson_propagator
#endif

#define GAMMAFIVE -1    // Some integer which is not a direction

// Flags for selecting M or Mdag
#define PLUS 1
#define MINUS -1

/*
* ROUTINES FOR SU(3) MATRIX OPERATIONS
*
* fundamental rep:
*
* void adjoint_f(matrix_f *a, matrix_f *b)
* file "adjoint_f.c"
* void scalar_add_diag_f(matrix_f *a, Real s)
* file "s_a_d_mat_f.c"
* void c_scalar_add_diag_f(matrix_f *a, complex *f)
* file "cs_a_d_mat_f.c"
* void make_anti_hermitian(matrix_f *m3,  anti_hermitmat *ah3)
* file "make_ahmat.c"
* void random_anti_hermitian(anti_hermitmat *mat_antihermit, double_prn *prn_pt)
* (prn_pt passed through to myrand())
* file "rand_ahmat.c"
* void uncompress_anti_hermitian(anti_hermitmat *mat_anti, matrix_f *mat)
* file "uncmp_ahmat.c"
* void compress_anti_hermitian(matrix_f *mat, anti_hermitmat *mat_anti)
* file "cmp_ahmat.c"
*
* fermion rep:
*
* complex det(matrix *a)
* file "det.c"
* void adjoint(matrix *a, matrix *b)
* file "adjoint.c"
*
* ROUTINES FOR IRREP VECTOR OPERATIONS (NCOL COMPONENT COMPLEX)
*
* Real rdot(vector *a, vector *b)
* file "rdot.c"
  *
* void mult_mat_vec(matrix *a, vector *b, vector *c)
*  C  <-  A*B
* file "m_matvec.c"
  *
* void mult_adj_mat_vec(matrix *a, vector *b, vector *c)
* file "m_amatvec.c"
*
* void sub_four_vecs(vector *a, vector *b1, vector *b2,
*   vector *b3, vector *b4)
* file "sub4vecs.c"
*
* ROUTINES MIXING SU(2) and SU(3)
*
* void left_su2_hit_n(su2_matrix *u, int p, int q, matrix *link)
*       file "l_su2_hit_n.c"
* void left_su2_hit_n_f(su2_matrix *u, int p, int q, matrix_f *link)
*       file "l_su2_hit_n_f.c"
* void right_su2_hit_a(su2_matrix *u, int p, int q, matrix *link)
*       file "r_su2_hit_a.c"
* void right_su2_hit_a_f(su2_matrix *u, int p, int q, matrix_f *link)
*       file "r_su2_hit_a_f.c"
* void dump_su2(su2_matrix *u)
*       file "dump_su2.c"
* void mult_su2_mat_vec_elem_n(su2_matrix *u, complex *x0, complex *x1)
*       file "m_su2_mat_vec_n.c"
* void mult_su2_mat_vec_elem_a(su2_matrix *u, complex *x0, complex *x1);
*       file "m_su2_mat_vec_a.c"
*
* ROUTINES FOR WILSON VECTORS
*
* void mult_mat_hwvec(matrix *mat, half_wilson_vector *src,
* half_wilson_vector *dest);
* file m_mat_hwvec.c
*    dest <- mat*src
* void mult_adj_mat_hwvec matrix *mat,
* half_wilson_vector *src, half_wilson_vector *dest)
* file m_amat_hwvec.c
*    dest <- mat_adjoint*src
*
* void wp_shrink(wilson_vector *src, half_wilson_vector *dest,
* int dir, int sign);
* file "wp_shrink.c"
*    if(dir = [XYZT]UP) dest <- components of src along eigenvectors
*   of gamma_dir with eigenvalue +1
*    if(dir = [XYZT]DOWN) dest <- components of src along eigenvectors
*   of gamma_dir with eigenvalue -1
*    if(sign==MINUS)switch roles of +1 and -1
* void wp_shrink_4dir(wilson_vector *a,  half_wilson_vector *b1,
* half_wilson_vector *b2, half_wilson_vector *b3,
* half_wilson_vector *b4, int sign);
* file "wp_shrink4.c"
*   Shrink A in X, Y, Z, T directions respectively, results in B1-B4
* void wp_grow(half_wilson_vector *src, wilson_vector *dest,
* int dir, int sign);
* file "wp_grow.c"
*    if(dir = [XYZT]UP) dest <- components of src times eigenvectors
*   of gamma_dir with eigenvalue +1
*    if(dir = [XYZT]DOWN) dest <- components of src times eigenvectors
*   of gamma_dir with eigenvalue -1
*    if(sign==MINUS)switch roles of +1 and -1
*       Note: wp_shrink(+-dir) followed by wp_grow(+-dir) amounts to
*   multiplication by 1+-gamma_dir, or 1-+gamma_dir if sign=MINUS
* void wp_grow_add(half_wilson_vector *src, wilson_vector *dest,
* int dir, int sign);
* file "wp_grow_a.c"
*    wp_grow, and add result to previous contents of dest.
* void grow_add_four_wvecs(wilson_vector *a, half_wilson_vector *b1,
* half_wilson_vector *b2, half_wilson_vector *b3,
* half_wilson_vector *b4, int sign, int sum);
* file "grow4wvecs.c"
*         If sum==0
*   Grow b1-b4 in X, Y, Z, T directions respectively, sum of results in A
*         If sum==1
*   Grow b1-b4 in X, Y, Z, T directions respectively, add to current A
*
* void mult_by_gamma(wilson_vector *src, wilson_vector *dest, int dir);
* file mb_gamma.c
*    dest <- gamma[dir] * src,  dir=[XYZT]UP, GAMMAFIVE
* void mult_by_gamma_left(wilson_matrix *src,  wilson_matrix *dest, int dir);
* file mb_gamma_l.c
*    dest <- gamma[dir] * src,  dir=[XYZT]UP, GAMMAFIVE
*    acts on first index of matrix
* void mult_by_gamma_right(wilson_matrix *src,  wilson_matrix *dest, int dir);
* file mb_gamma_r.c
*    dest_ij <- gamma[dir]_ik * src_jk,  dir=[XYZT]UP, GAMMAFIVE
*    acts on second index of matrix
*
* void mult_swv_by_gamma_l(spin_wilson_vector *src,  spin_wilson_vector *dest, int dir);
* file mswvb_gamma_l.c
*    same as mult_by_gamma_left, but acts on a spin_wilson_vector
*
* void mult_swv_by_gamma_r(spin_wilson_vector *src,  spin_wilson_vector *dest, int dir);
* file mswvb_gamma_r.c
*    same as mult_by_gamma_right, but acts on a spin_wilson_vector
*
* void projector_w(wilson_vector *a, wilson_vector *b, matrix *c)
* sum over spins of outer product of A.d[s] and B.d[s]  - a three
*   by three complex matrix
* file "proj_w.c"
* void copy_wvec(wilson_vector *src, wilson_vector *dest);
* file copy_wvec.c
*    dest <- src
*/
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Subroutine definitions
// Fundamental vector operations
// In file clear_vec_f.c
void clear_vec_f(vector_f *v);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Irrep vector operations
// In file dump_vec.c
void dump_vec(vector *v);

// In file clear_vec.c
void clear_vec(vector *v);

// In file vec_copy.c
void vec_copy(vector *a, vector *b);

// In file msq_su3vec.c
Real magsq_su3vec(vector *a);
void magsq_su3vec_sum(vector *a, Real *c);

// In file addvec.c
void sum_vector(vector *b, vector *c);
void dif_vector(vector *b, vector *c);
void add_vector(vector *a, vector *b, vector *c);
void sub_vector(vector *a, vector *b, vector *c);

// In file s_m_vec.c
void scalar_mult_vector(vector *b, Real s, vector *c);
void scalar_mult_sum_vector(vector *b, Real s, vector *c);
void scalar_mult_dif_vector(vector *b, Real s, vector *c);
void scalar_mult_add_vector(vector *a, vector *b, Real s, vector *c);

// In file cs_m_a_vec.c
void c_scalar_mult_add_su3vec(vector *a, vector *b, complex *s, vector *c);
void c_scalar_mult_sum_su3vec(vector *b, complex *s, vector *c);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Anti-hermitian matrix routines
void make_anti_hermitian(matrix_f *m, anti_hermitmat *ah);
void random_anti_hermitian(anti_hermitmat *ah, double_prn *prn_pt);
void uncompress_anti_hermitian(anti_hermitmat *ah, matrix_f *m);
void compress_anti_hermitian(matrix_f *m, anti_hermitmat *ah);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Fundamental matrix operations
// In file dump_mat_f.c
void dump_mat_f(matrix_f *m);

// In file clear_mat_f.c
void clear_mat_f(matrix_f *m);

// In file trace_f.c
complex trace_f(matrix_f *a);
void trace_sum_f(matrix_f *a, complex *c);

// In file tr_prod_f.c
Real realtrace_f(matrix_f *a, matrix_f *b);
Real realtrace_nn_f(matrix_f *a, matrix_f *b);
void realtrace_sum_f(matrix_f *a, matrix_f *b, Real *c);
complex complextrace_f(matrix_f *a, matrix_f *b);

// In file mat_copy_f.c
void mat_copy_f(matrix_f *a, matrix_f *b);

// In file addmat_f.c
void sum_mat_f(matrix_f *b, matrix_f *c);
void dif_mat_f(matrix_f *b, matrix_f *c);
void add_mat_f(matrix_f *a, matrix_f *b, matrix_f *c);
void sub_mat_f(matrix_f *a, matrix_f *b, matrix_f *c);

// In file s_m_mat_f.c
void scalar_mult_mat_f(matrix_f *src, Real scalar, matrix_f *dest);
void scalar_mult_sum_mat_f(matrix_f *b, Real s, matrix_f *c);
void scalar_mult_dif_mat_f(matrix_f *b, Real s, matrix_f *c);
void scalar_mult_add_mat_f(matrix_f *a, matrix_f *b, Real s, matrix_f *c);

// In file cs_m_mat_f.c
void c_scalar_mult_mat_f(matrix_f *b, complex *s, matrix_f *c);
void c_scalar_mult_sum_mat_f(matrix_f *b, complex *s, matrix_f *c);
void c_scalar_mult_add_mat_f(matrix_f *a, matrix_f *b, complex *s,
                             matrix_f *c);

// In file m_mat_nn_f.c
void mult_nn_sum_f(matrix_f *a, matrix_f *b, matrix_f *c);
void mult_nn_f(matrix_f *a, matrix_f *b, matrix_f *c);

// In file m_mat_na_f.c
void mult_na_sum_f(matrix_f *a, matrix_f *b, matrix_f *c);
void mult_na_f(matrix_f *a, matrix_f *b, matrix_f *c);

// In file m_mat_an_f.c
void mult_an_sum_f(matrix_f *a, matrix_f *b, matrix_f *c);
void mult_an_f(matrix_f *a, matrix_f *b, matrix_f *c);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Irrep matrix operations
// In file dump_mat.c
void dump_mat(matrix *m);

// In file clear_mat.c
void clear_mat(matrix *m);

// In file trace.c
complex trace(matrix *a);

// In file tr_prod.c
Real realtrace(matrix *a, matrix *b);
void realtrace_sum(matrix *a, matrix *b, Real *c);
complex complextrace(matrix *a, matrix *b);

// In file mat_copy.c
void mat_copy(matrix *a, matrix *b);

// In file addmat.c
void sum_mat(matrix *b, matrix *c);
void dif_mat(matrix *b, matrix *c);
void add_mat(matrix *a, matrix *b, matrix *c);
void sub_mat(matrix *a, matrix *b, matrix *c);

// In file s_m_mat.c
void scalar_mult_mat(matrix *src, Real scalar, matrix *dest);
void scalar_mult_sum_mat(matrix *b, Real s, matrix *c);
void scalar_mult_add_mat(matrix *a, matrix *b, Real s, matrix *c);

// In file m_mat_nn.c
void mult_nn_dif(matrix *a, matrix *b, matrix *c);
void mult_nn(matrix *a, matrix *b, matrix *c);

// In file m_mat_na.c
void mult_na_sum(matrix *a, matrix *b, matrix *c);
void mult_na_dif(matrix *a, matrix *b, matrix *c);
void mult_na(matrix *a, matrix *b, matrix *c);

// In file m_mat_an.c
void mult_an(matrix *a, matrix *b, matrix *c);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Wilson vector operations
// In file dump_wvec.c
void dump_wvec(wilson_vector *src);

// In file clear_wvec.c
void clear_wvec(wilson_vector *dest);

// In file msq_wvec.c
Real magsq_wvec(wilson_vector *a);
void magsq_wvec_sum(wilson_vector *a, Real *c);

// In file wvec_dot.c
Real wvec_rdot(wilson_vector *a, wilson_vector *b);
complex wvec_dot(wilson_vector *a, wilson_vector *b);
void wvec_rdot_sum(wilson_vector *a, wilson_vector *b, Real *c);
void wvec_dot_sum(wilson_vector *a, wilson_vector *b, complex *c);

// In file addwvec.c
void sum_wvec(wilson_vector *b, wilson_vector *c);
void dif_wvec(wilson_vector *b, wilson_vector *c);
void add_wvec(wilson_vector *a, wilson_vector *b, wilson_vector *c);
void sub_wvec(wilson_vector *a, wilson_vector *b, wilson_vector *c);

// In file s_m_wvec.c
void scalar_mult_wvec(wilson_vector *b, Real s, wilson_vector *c);
void scalar_mult_sum_wvec(wilson_vector *b, Real s, wilson_vector *c);
void scalar_mult_dif_wvec(wilson_vector *b, Real s, wilson_vector *c);
void scalar_mult_add_wvec(wilson_vector *a, wilson_vector *b, Real s,
                          wilson_vector *c);

// In file cs_m_wvec.c
void c_scalar_mult_sum_wvec(wilson_vector *b, complex *s, wilson_vector *c);
void c_scalar_mult_add_wvec(wilson_vector *a, wilson_vector *b, complex *s,
                            wilson_vector *c);

// In file m_mat_wvec.c
void mult_mat_wvec(matrix *a, wilson_vector *b, wilson_vector *c);
void mult_adj_mat_wvec(matrix *a, wilson_vector *b, wilson_vector *c);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Routines mixing SU(2) and U(N)
void left_su2_hit_n(su2_matrix *u, int p, int q, matrix *link);
void right_su2_hit_a(su2_matrix *u, int p, int q, matrix *link);
void left_su2_hit_n_f(su2_matrix *u, int p, int q, matrix_f *link);
void right_su2_hit_a_f(su2_matrix *u, int p, int q, matrix_f *link);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Miscellaneous routines
// In file gaussrand.c
Real gaussian_rand_no(double_prn *prn_pt);

// In file byterevn.c
#include "../include/int32type.h"
void byterevn(int32type w[], int n);
void byterevn64(int32type w[], int n);
// -----------------------------------------------------------------


void adjoint_f(matrix_f *a, matrix_f *b);

complex det(matrix *a);
void adjoint(matrix *a, matrix *b);

void scalar_add_diag_f(matrix_f *a, Real s);
void c_scalar_add_diag_f(matrix_f *a, complex *s);

void dump_su2(su2_matrix *u);
void mult_su2_mat_vec_elem_n(su2_matrix *u, complex *x0, complex *x1);
void mult_su2_mat_vec_elem_a(su2_matrix *u, complex *x0, complex *x1);

void wp_shrink(wilson_vector *src, half_wilson_vector *dest,
               int dir, int sign);

// In file wp_grow.c
void wp_grow(half_wilson_vector *src, wilson_vector *dest,
             int dir, int sign);
void wp_grow_sum(half_wilson_vector *src, wilson_vector *dest,
                 int dir, int sign);

// In file grow4wvecs.c
void grow_add_four_wvecs(wilson_vector *a, half_wilson_vector *b1,
                         half_wilson_vector *b2, half_wilson_vector *b3,
                         half_wilson_vector *b4, int sign);
void grow_sum_four_wvecs(wilson_vector *a, half_wilson_vector *b1,
                         half_wilson_vector *b2, half_wilson_vector *b3,
                         half_wilson_vector *b4, int sign);

void mult_by_gamma(wilson_vector *src, wilson_vector *dest, int dir);
void mult_by_gamma_left(wilson_matrix *src,  wilson_matrix *dest, int dir);
void mult_by_gamma_right(wilson_matrix *src,  wilson_matrix *dest, int dir);
void mult_swv_by_gamma_l(spin_wilson_vector *src, spin_wilson_vector *dest,
                         int dir);
void mult_swv_by_gamma_r(spin_wilson_vector *src, spin_wilson_vector *dest,
                         int dir);
void projector_w(wilson_vector *a, wilson_vector *b, matrix *c);
void copy_wvec(wilson_vector *src, wilson_vector *dest);

void mult_mat_vec(matrix *a, vector *b, vector *c);

void mult_adj_mat_vec(matrix *a, vector *b, vector *c);

void mult_adj_mat_hwvec(matrix *mat, half_wilson_vector *src,
                        half_wilson_vector *dest);

void mult_mat_hwvec(matrix *mat, half_wilson_vector *src,
                    half_wilson_vector *dest);

void sub_four_vecs(vector *a, vector *b1, vector *b2,
                   vector *b3, vector *b4);

void wp_shrink_4dir(wilson_vector *a,  half_wilson_vector *b1,
                    half_wilson_vector *b2, half_wilson_vector *b3,
                    half_wilson_vector *b4, int sign);

#endif
// -----------------------------------------------------------------
