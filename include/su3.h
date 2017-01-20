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
// FREP can take the values: fundamental, symmetric2, antisymmetric2
#define NCOL 3
#define DIMF 3
#define FREP fundamental

// Only have anti_hermitmat and explicit loops in libraries for N <= 4
#if NCOL > 4
  #error "NCOL > 4 not yet implemented!"
#endif

#if (NCOL != 3 || DIMF != 3)
  #ifdef FAST
    #error "FAST only works if NCOL = DIMF = 3"
  #endif
#endif

typedef struct { fcomplex e[NCOL][NCOL]; } fsu3_matrix_f;
typedef struct { fcomplex e[DIMF][DIMF]; } fsu3_matrix;
typedef struct { fcomplex c[NCOL]; } fsu3_vector_f;
typedef struct { fcomplex c[DIMF]; } fsu3_vector;
typedef struct {
#if NCOL == 2
  fcomplex m01;
  float m00im, m11im;
#endif
#if NCOL == 3
  fcomplex m01, m02, m12;
  float m00im, m11im, m22im;
  float space;
#endif
#if NCOL == 4
  fcomplex m01, m02, m03, m12, m13, m23;
  float m00im, m11im, m22im, m33im;
#endif
} fanti_hermitmat;

typedef struct { dcomplex e[NCOL][NCOL]; } dsu3_matrix_f;
typedef struct { dcomplex e[DIMF][DIMF]; } dsu3_matrix;
typedef struct { dcomplex c[NCOL]; } dsu3_vector_f;
typedef struct { dcomplex c[DIMF]; } dsu3_vector;
typedef struct {
#if NCOL == 2
  dcomplex m01;
  double m00im, m11im;
#endif
#if NCOL == 3
  dcomplex m01, m02, m12;
  double m00im, m11im, m22im;
  double space;
#endif
#if NCOL == 4
  dcomplex m01, m02, m03, m12, m13, m23;
  double m00im, m11im, m22im, m33im;
#endif
} danti_hermitmat;

#if PRECISION == 1
#define su3_matrix_f   fsu3_matrix_f
#define su3_matrix     fsu3_matrix
#define su3_vector_f   fsu3_vector_f
#define su3_vector     fsu3_vector
#define anti_hermitmat fanti_hermitmat
#else
#define su3_matrix_f   dsu3_matrix_f
#define su3_matrix     dsu3_matrix
#define su3_vector_f   dsu3_vector_f
#define su3_vector     dsu3_vector
#define anti_hermitmat danti_hermitmat
#endif

/* SU(2) */
typedef struct { complex e[2][2]; } su2_matrix;

/* Wilson vectors */
/* e.g.                */
/* wilson_propagator prop;                           */
/* prop.c[ci].d[si].d[sf].c[cf]                      */
/* ----------------------->    complex               */
/* ----------------->          su3_vector            */
/* ----------->                wilson_vector         */
/* ----->                      spin_wilson_vector    */
/* e.g.                */
/* wilson_matrix matr;                               */
/* matr.d[si].c[ci].d[sf].c[cf]                      */
/* ----------------------->    complex               */
/* ----------------->          su3_vector            */
/* ----------->                wilson_vector         */
/* ----->                      color_wilson_vector   */

/* Object with two Dirac and two color indices. A given element
   of a "wilson_propagator" is accessed by
   object.c[color1].d[spin1].d[spin2].c[color2].real , etc.
   As alway, "d" denotes a Dirac index and "c" a color index.
   "1" refers to the source, "2" to the sink.

   Color indices here always run to DIMF, not NCOL
*/

typedef struct { fsu3_vector d[4]; } fwilson_vector;
typedef struct { fsu3_vector h[2]; } fhalf_wilson_vector;
typedef struct { fwilson_vector c[DIMF]; } fcolor_wilson_vector;
typedef struct { fwilson_vector d[4]; } fspin_wilson_vector;
typedef struct { fcolor_wilson_vector d[4]; } fwilson_matrix;
typedef struct { fspin_wilson_vector c[DIMF]; } fwilson_propagator;

typedef struct { dsu3_vector d[4]; } dwilson_vector;
typedef struct { dsu3_vector h[2]; } dhalf_wilson_vector;
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
* void mult_su3_an_f(su3_matrix_f *a, su3_matrix_f *b, su3_matrix_f *c)
* file "m_mat_an_f.c"
* void mult_su3_na_f(su3_matrix_f *a, su3_matrix_f *b, su3_matrix_f *c)
* file "m_mat_na_f.c"
* void mult_su3_nn_f(su3_matrix_f *a, su3_matrix_f *b, su3_matrix_f *c)
* file "m_mat_nn_f.c"
* void add_su3_matrix_f(su3_matrix_f *a, su3_matrix_f *b, su3_matrix_f *c);
*       file "addmat_f.c"
* Real realtrace_su3_f(su3_matrix_f *a, su3_matrix_f *b)
* file "realtr_f.c"
* void scalar_mult_add_su3_matrix_f(su3_matrix_f *a, su3_matrix_f *b,
*       Real s, su3_matrix_f *c)
* file s_m_a_mat_f.c
* void scalar_mult_sub_su3_matrix_f(su3_matrix_f *src1, su3_matrix_f *src2,
  Real scalar, su3_matrix_f *dest);
* file "s_m_s_mat_f.c"
* void su3_adjoint_f(su3_matrix_f *a, su3_matrix_f *b)
* file "su3_adjoint_f.c"
* void su3mat_copy_f(su3_matrix_f *a, su3_matrix_f *b)
* file "su3mat_copy_f.c"
* void clearvec_f(su3_vector_f *v)
*       file "clearvec_f.c"
* complex trace_su3_f(su3_matrix_f *a)
* file "trace_su3_f.c"
* void sub_su3_matrix_f(su3_matrix_f *a, su3_matrix_f *b, su3_matrix_f *c)
* file "submat_f.c"
* void scalar_mult_su3_matrix_f(su3_matrix_f *a, Real s, su3_matrix_f *b)
* file "s_m_mat_f.c"
* void scalar_add_diag_su3_f(su3_matrix_f *a, Real s)
* file "s_a_d_mat_f.c"
* void c_scalar_add_diag_su3_f(su3_matrix_f *a, complex *f)
* file "cs_a_d_mat_f.c"
* complex complextrace_su3_f(su3_matrix_f *a, su3_matrix_f *b)
*       file "complextr_f.c"
* void c_scalar_mult_su3mat_f(su3_matrix_f *b, complex *s, su3_matrix_f *c)
*       file "cs_m_mat_f.c"
* void c_scalar_mult_add_su3mat_f(su3_matrix_f *m1, su3_matrix_f *m2,
*        complex *phase, su3_matrix_f *m3)
*       file "cs_m_a_mat_f.c"
* void make_anti_hermitian(su3_matrix_f *m3,  anti_hermitmat *ah3)
* file "make_ahmat.c"
* void random_anti_hermitian(anti_hermitmat *mat_antihermit, double_prn *prn_pt)
* (prn_pt passed through to myrand())
* file "rand_ahmat.c"
* void uncompress_anti_hermitian(anti_hermitmat *mat_anti, su3_matrix_f *mat)
* file "uncmp_ahmat.c"
* void compress_anti_hermitian(su3_matrix_f *mat, anti_hermitmat *mat_anti)
* file "cmp_ahmat.c"
* void dumpmat_f(su3_matrix_f *m)
*       file "dumpmat_f.c"
*
* fermion rep:
*
* void mult_su3_nn(su3_matrix *a, su3_matrix *b, su3_matrix *c )
* matrix multiply, no adjoints
* files "m_mat_nn.c"
* void mult_su3_na(su3_matrix *a, su3_matrix *b, su3_matrix *c)
* matrix multiply, second matrix is adjoint
* files "m_mat_na.c"
* void mult_su3_an(su3_matrix *a, su3_matrix *b, su3_matrix *c)
* matrix multiply, first matrix is adjoint
* files "m_mat_an.c"
* Real realtrace_su3(su3_matrix *a, su3_matrix *b)
* (Re(Tr(A_adjoint*B)))
* file "realtr.c"
* complex trace_su3(su3_matrix *a)
* file "trace_su3.c"
* complex complextrace_su3(su3_matrix *a, su3_matrix *b)
* (Tr(A_adjoint*B))
* file "complextr.c"
* complex det_su3(su3_matrix *a)
* file "det_su3.c"
* void add_su3_matrix(su3_matrix *a, su3_matrix *b, su3_matrix *c)
* file "addmat.c"
* void sub_su3_matrix(su3_matrix *a, su3_matrix *b, su3_matrix *c)
* file "submat.c"
* void scalar_mult_su3_matrix(su3_matrix *a, Real s, su3_matrix *b)
* file "s_m_mat.c"
* void scalar_mult_add_su3_matrix(su3_matrix *a, su3_matrix *b,
* Real s, su3_matrix *c)
* file "s_m_a_mat.c"
* void su3_adjoint(su3_matrix *a, su3_matrix *b)
* file "su3_adjoint.c"
* void clear_su3mat(su3_matrix *dest);
*       file clear_mat.c
* void su3mat_copy(su3_matrix *a, su3_matrix *b)
* file "su3mat_copy.c"
* void dumpmat(su3_matrix *m)
*       file "dumpmat.c"
*
* ROUTINES FOR IRREP VECTOR OPERATIONS (NCOL COMPONENT COMPLEX)
*
* Real su3_rdot(su3_vector *a, su3_vector *b)
* file "su3_rdot.c"
  *
* Real magsq_su3vec(su3_vector *a)
* file "msq_su3vec.c"
  *
* void su3vec_copy(su3_vector *a, su3_vector *b)
* file "su3vec_copy.c"
*
* void mult_su3_mat_vec(su3_matrix *a, su3_vector *b, su3_vector *c)
*  C  <-  A*B
* file "m_matvec.c"
  *
* void mult_adj_su3_mat_vec(su3_matrix *a, su3_vector *b, su3_vector *c)
* file "m_amatvec.c"
*
* void add_su3_vector(su3_vector *a, su3_vector *b, su3_vector *c)
* file "addvec.c"
* void sub_su3_vector(su3_vector *a, su3_vector *b, su3_vector *c)
* file "subvec.c"
* void sub_four_su3_vecs(su3_vector *a, su3_vector *b1, su3_vector *b2,
*   su3_vector *b3, su3_vector *b4)
* file "sub4vecs.c"
*
* void scalar_mult_su3_vector(su3_vector *a, Real s, su3_vector *c)
* file "s_m_vec.c"
* void scalar_mult_add_su3_vector(su3_vector *a, su3_vector *b, Real s,
* su3_vector *c)
* file "s_m_a_vec.c"
* void scalar_mult_sum_su3_vector(su3_vector *a, su3_vector *b, Real s)
* file "s_m_sum_vec.c"
* void dumpvec(su3_vector *v)
*       file "dumpvec.c"
* void clearvec(su3_vector *v)
*       file "clearvec.c"
*
* ROUTINES MIXING SU(2) and SU(3)
*
* void left_su2_hit_n(su2_matrix *u, int p, int q, su3_matrix *link)
*       file "l_su2_hit_n.c"
* void left_su2_hit_n_f(su2_matrix *u, int p, int q, su3_matrix_f *link)
*       file "l_su2_hit_n_f.c"
* void right_su2_hit_a(su2_matrix *u, int p, int q, su3_matrix *link)
*       file "r_su2_hit_a.c"
* void right_su2_hit_a_f(su2_matrix *u, int p, int q, su3_matrix_f *link)
*       file "r_su2_hit_a_f.c"
* void dumpsu2(su2_matrix *u)
*       file "dumpsu2.c"
* void mult_su2_mat_vec_elem_n(su2_matrix *u, complex *x0, complex *x1)
*       file "m_su2_mat_vec_n.c"
* void mult_su2_mat_vec_elem_a(su2_matrix *u, complex *x0, complex *x1);
*       file "m_su2_mat_vec_a.c"
*
* ROUTINES FOR WILSON VECTORS
*
* void mult_mat_wilson_vec(su3_matrix *mat, wilson_vector *src,
* wilson_vector *dest);
* file m_mat_wvec.c
*    dest <- mat*src
* void mult_su3_mat_hwvec(su3_matrix *mat, half_wilson_vector *src,
* half_wilson_vector *dest);
* file m_mat_hwvec.c
*    dest <- mat*src
* void mult_adj_mat_wilson_vec(su3_matrix *mat, wilson_vector *src,
* wilson_vector *dest)
* file m_amat_wvec.c
*    dest <- mat_adjoint*src
* void mult_adj_su3_mat_hwvec su3_matrix *mat,
* half_wilson_vector *src, half_wilson_vector *dest)
* file m_amat_hwvec.c
*    dest <- mat_adjoint*src
*
* void add_wilson_vector(wilson_vector *src1, wilson_vector *src2,
* wilson_vector *dest);
* file add_wvec.c
*    dest <- src1+src2
* void sub_wilson_vector(wilson_vector *src1, wilson_vector *src2,
*       wilson_vector *dest);
* file sub_wvec.c
*    dest <- src1-src2
*
* void scalar_mult_wvec  wilson_vector *src, Real s, wilson_vector *dest)
* file s_m_wvec.c
*    dest <- s*src
* Real magsq_wvec(wilson_vector *src);
* file msq_wvec.c
*    s <- squared magnitude of src
* complex wvec_dot(wilson_vector *src1, wilson_vector *src2);
* file wvec_dot.c
*    c <- dot product of src1 and src2
* Real wvec_rdot(wilson_vector *a, wilson_vector *b)
* wilson_vector *a,*b;
* file "wvec_rdot.c"
*    r <- real part of dot product of src1 and src2
* void scalar_mult_add_wvec(wilson_vector *src1,*src2, Real s,
*           wilson_vector *dest)
* file s_m_a_wvec.c
*    dest <- src1 + s*src2
* void c_scalar_mult_add_wvec(wilson_vector *v1, wilson_vector *v2,
* complex phase, wilson_vector *v3)
*       file "cs_m_a_wvec.c"
* differs from previous one: value of scalar, not address is arg.
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
* void su3_projector_w(wilson_vector *a, wilson_vector *b, su3_matrix *c)
* sum over spins of outer product of A.d[s] and B.d[s]  - a three
*   by three complex matrix
* file "su3_proj_w.c"
* void clear_wvec(wilson_vector *dest);
* file clear_wvec.c
*    dest <- 0.0
* void copy_wvec(wilson_vector *src, wilson_vector *dest);
* file copy_wvec.c
*    dest <- src
* void dump_wvec(wilson_vector *src);
* file dump_wvec.c
*    print out a wilson vector
*
* MISCELLANEOUS ROUTINES
*
* Real gaussian_rand_no(double_prn *prn_pt)
* file "gaussrand.c"
* void byterevn(int32type w[], int n)
* void byterevn64(int32type w[], int n)
*
*/

Real realtrace_su3_f(su3_matrix_f *a, su3_matrix_f *b);

void su3_adjoint_f(su3_matrix_f *a, su3_matrix_f *b);
void su3mat_copy_f(su3_matrix_f *a, su3_matrix_f *b);
complex trace_su3_f(su3_matrix_f *a);
void clear_su3mat_f(su3_matrix_f *dest);
void clearvec_f(su3_vector_f *v);
Real realtrace_su3(su3_matrix *a, su3_matrix *b);
complex trace_su3(su3_matrix *a);
complex complextrace_su3(su3_matrix *a, su3_matrix *b);
complex det_su3(su3_matrix *a);
void scalar_mult_su3_matrix(su3_matrix *src, Real scalar, su3_matrix *dest);
void scalar_mult_su3_matrix_f(su3_matrix_f *src, Real scalar, su3_matrix_f *dest);
void su3_adjoint(su3_matrix *a, su3_matrix *b);
void scalar_add_diag_su3_f(su3_matrix_f *a, Real s);
void c_scalar_add_diag_su3_f(su3_matrix_f *a, complex *f);
complex complextrace_su3_f(su3_matrix_f *a, su3_matrix_f *b);

// In file cs_m_mat_f.c
void c_scalar_mult_su3mat_f(su3_matrix_f *b, complex *s, su3_matrix_f *c);
void c_scalar_mult_sum_su3mat_f(su3_matrix_f *b, complex *s, su3_matrix_f *c);
void c_scalar_mult_add_su3mat_f(su3_matrix_f *a, su3_matrix_f *b, complex *s,
                                su3_matrix_f *c);

void make_anti_hermitian(su3_matrix_f *m3, anti_hermitmat *ah3);
void random_anti_hermitian(anti_hermitmat *mat_antihermit, double_prn *prn_pt);
void uncompress_anti_hermitian(anti_hermitmat *mat_anti, su3_matrix_f *mat);
void compress_anti_hermitian(su3_matrix_f *mat, anti_hermitmat *mat_anti);
void clear_su3mat(su3_matrix *dest);
void su3mat_copy(su3_matrix *a, su3_matrix *b);
void dumpmat(su3_matrix *m);
void dumpmat_f(su3_matrix_f *m);

void su3vec_copy(su3_vector *a, su3_vector *b);
void dumpvec(su3_vector *v);
void clearvec(su3_vector *v);

// In file addvec.c
void sum_su3_vector(su3_vector *b, su3_vector *c);
void dif_su3_vector(su3_vector *b, su3_vector *c);
void add_su3_vector(su3_vector *a, su3_vector *b, su3_vector *c);
void sub_su3_vector(su3_vector *a, su3_vector *b, su3_vector *c);

// In file s_m_vec.c
void scalar_mult_su3_vector(su3_vector *b, Real s, su3_vector *c);
void scalar_mult_sum_su3_vector(su3_vector *b, Real s, su3_vector *c);
void scalar_mult_dif_su3_vector(su3_vector *b, Real s, su3_vector *c);
void scalar_mult_add_su3_vector(su3_vector *a, su3_vector *b, Real s,
                                su3_vector *c);

// In file cs_m_a_vec.c
void c_scalar_mult_add_su3vec(su3_vector *a, su3_vector *b, complex *s,
                              su3_vector *c);
void c_scalar_mult_sum_su3vec(su3_vector *b, complex *s, su3_vector *c);

// In file wvec_dot.c
Real wvec_rdot(wilson_vector *a, wilson_vector *b);
complex wvec_dot(wilson_vector *a, wilson_vector *b);
void wvec_rdot_sum(wilson_vector *a, wilson_vector *b, Real *c);
void wvec_dot_sum(wilson_vector *a, wilson_vector *b, complex *c);

// In file addwvec.c
void sum_wvec(wilson_vector *b, wilson_vector *c);
void dif_wvec(wilson_vector *b, wilson_vector *c);
void add_wilson_vector(wilson_vector *a, wilson_vector *b, wilson_vector *c);
void sub_wilson_vector(wilson_vector *a, wilson_vector *b, wilson_vector *c);

// In file s_m_wvec.c
void scalar_mult_wvec(wilson_vector *b, Real s, wilson_vector *c);
void scalar_mult_sum_wvec(wilson_vector *b, Real s, wilson_vector *c);
void scalar_mult_dif_wvec(wilson_vector *b, Real s, wilson_vector *c);
void scalar_mult_add_wvec(wilson_vector *a, wilson_vector *b, Real s,
                          wilson_vector *c);

// In file cs_m_a_wvec.c
void c_scalar_mult_sum_wvec(wilson_vector *b, complex *s, wilson_vector *c);
void c_scalar_mult_add_wvec(wilson_vector *a, wilson_vector *b, complex *s,
                            wilson_vector *c);

void left_su2_hit_n(su2_matrix *u, int p, int q, su3_matrix *link);
void right_su2_hit_a(su2_matrix *u, int p, int q, su3_matrix *link);
void left_su2_hit_n_f(su2_matrix *u, int p, int q, su3_matrix_f *link);
void right_su2_hit_a_f(su2_matrix *u, int p, int q, su3_matrix_f *link);
void dumpsu2(su2_matrix *u);
void mult_su2_mat_vec_elem_n(su2_matrix *u, complex *x0, complex *x1);
void mult_su2_mat_vec_elem_a(su2_matrix *u, complex *x0, complex *x1);

void mult_mat_wilson_vec(su3_matrix *mat, wilson_vector *src,
                         wilson_vector *dest);
void mult_adj_mat_wilson_vec(su3_matrix *mat, wilson_vector *src,
                             wilson_vector *dest);

// In file msq_wvec.c
Real magsq_wvec(wilson_vector *src);
void magsq_wvec_sum(wilson_vector *src, Real *c);

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
void su3_projector_w(wilson_vector *a, wilson_vector *b, su3_matrix *c);
void clear_wvec(wilson_vector *dest);
void copy_wvec(wilson_vector *src, wilson_vector *dest);
void dump_wvec(wilson_vector *src);

Real gaussian_rand_no(double_prn *prn_pt);
#include "../include/int32type.h"
void byterevn(int32type w[], int n);
void byterevn64(int32type w[], int n);

// In file addmat.c
void sum_su3_matrix(su3_matrix *b, su3_matrix *c);
void dif_su3_matrix(su3_matrix *b, su3_matrix *c);
void add_su3_matrix(su3_matrix *a, su3_matrix *b, su3_matrix *c);
void sub_su3_matrix(su3_matrix *a, su3_matrix *b, su3_matrix *c);

// In file addmat_f.c
void sum_su3_matrix_f(su3_matrix_f *b, su3_matrix_f *c);
void dif_su3_matrix_f(su3_matrix_f *b, su3_matrix_f *c);
void add_su3_matrix_f(su3_matrix_f *a, su3_matrix_f *b, su3_matrix_f *c);
void sub_su3_matrix_f(su3_matrix_f *a, su3_matrix_f *b, su3_matrix_f *c);

// In file msq_su3vec.c
Real magsq_su3vec(su3_vector *a);
void magsq_su3vec_sum(su3_vector *a, Real *c);

// In file m_mat_nn.c
void mult_su3_nn_dif(su3_matrix *a, su3_matrix *b, su3_matrix *c);
void mult_su3_nn(su3_matrix *a, su3_matrix *b, su3_matrix *c);

// In file m_mat_na.c
void mult_su3_na_sum(su3_matrix *a, su3_matrix *b, su3_matrix *c);
void mult_su3_na_dif(su3_matrix *a, su3_matrix *b, su3_matrix *c);
void mult_su3_na(su3_matrix *a, su3_matrix *b, su3_matrix *c);

// In file m_mat_an.c
void mult_su3_an(su3_matrix *a, su3_matrix *b, su3_matrix *c);

// In file m_mat_nn_f.c
void mult_su3_nn_f(su3_matrix_f *a, su3_matrix_f *b, su3_matrix_f *c);

// In file m_mat_na_f.c
void mult_su3_na_sum_f(su3_matrix_f *a, su3_matrix_f *b, su3_matrix_f *c);
void mult_su3_na_f(su3_matrix_f *a, su3_matrix_f *b, su3_matrix_f *c);

// In file m_mat_an_f.c
void mult_su3_an_sum_f(su3_matrix_f *a, su3_matrix_f *b, su3_matrix_f *c);
void mult_su3_an_f(su3_matrix_f *a, su3_matrix_f *b, su3_matrix_f *c);

void mult_su3_mat_vec(su3_matrix *a, su3_vector *b, su3_vector *c);

void mult_adj_su3_mat_vec(su3_matrix *a, su3_vector *b, su3_vector *c);

void mult_adj_su3_mat_hwvec(su3_matrix *mat, half_wilson_vector *src,
                            half_wilson_vector *dest);

void mult_su3_mat_hwvec(su3_matrix *mat, half_wilson_vector *src,
                        half_wilson_vector *dest);

// In file s_m_a_mat.c
void scalar_mult_sum_su3_matrix(su3_matrix *b, Real s, su3_matrix *c);
void scalar_mult_add_su3_matrix(su3_matrix *a, su3_matrix *b, Real s,
                                su3_matrix *c);

// In file s_m_a_mat_f.c
void scalar_mult_sum_su3_matrix_f(su3_matrix_f *b, Real s, su3_matrix_f *c);
void scalar_mult_add_su3_matrix_f(su3_matrix_f *a, su3_matrix_f *b, Real s,
                                  su3_matrix_f *c);

// In file s_m_s_mat_f.c
void scalar_mult_dif_su3_matrix_f(su3_matrix_f *b, Real s, su3_matrix_f *c);
void scalar_mult_sub_su3_matrix_f(su3_matrix_f *a, su3_matrix_f *b, Real s,
                                  su3_matrix_f *c);

void sub_four_su3_vecs(su3_vector *a, su3_vector *b1, su3_vector *b2,
  su3_vector *b3, su3_vector *b4);

void sub_su3_vector(su3_vector *a, su3_vector *b, su3_vector *c);

void wp_shrink_4dir(wilson_vector *a,  half_wilson_vector *b1,
  half_wilson_vector *b2, half_wilson_vector *b3,
  half_wilson_vector *b4, int sign);

#endif
// -----------------------------------------------------------------
