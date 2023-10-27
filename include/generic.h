// -----------------------------------------------------------------
// Macros and declarations for miscellaneous generic routines
#ifndef _GENERIC_H
#define _GENERIC_H

// Other generic directory declarations are elsewhere:
//   See comdefs.h for communications
//   See io_lat.h for I/O
#include <stdio.h>
#include "../include/int32type.h"
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/macros.h"
#include "../include/random.h"
#include "../include/file_types.h"
#include "../include/io_lat.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// check_unitarity.c
Real check_unit(matrix *c);
Real check_unitarity();

// nersc_cksum.c
void linktrsum(double_complex *linktr);

// plaquette.c
void plaquette(double *ss_plaq, double *st_plaq);

// field_strength.c
// link_src is offset for matrix link[4] in site struct
// field_dest is offset for matrix fieldstrength[6] in site struct
void make_field_strength(field_offset link_src, field_offset field_dest);

// gaugefix.c
void gaugefix(int gauge_dir, Real relax_boost, int max_gauge_iter,
              Real gfix_tol, field_offset diffmat, field_offset sumvec);

/* io_helpers.c */
gauge_file *save_lattice(int flag, char *filename, char *stringLFN );
gauge_file *reload_lattice(int flag, char *filename);
int ask_starting_lattice(FILE *fp, int prompt, int *flag, char *filename );
int ask_ending_lattice(FILE *fp, int prompt, int *flag, char *filename );
int ask_ildg_LFN(FILE *fp, int prompt, int flag, char *stringLFN);
void coldlat();
void funnylat();
int get_f(FILE *fp, int prompt, char *variable_name_string, Real *value );
int get_i(FILE *fp, int prompt, char *variable_name_string, int *value );
int get_vi(FILE *fp, int prompt, char *variable_name_string,
      int *value, int nvalues );
int get_vf(FILE *fp, int prompt, char *variable_name_string,
      Real *value, int nvalues );
int get_s(FILE *fp, int prompt, char *variable_name_string, char *value );
int get_prompt(FILE *fp, int *value );

/* layout_hyper_prime.c */
int io_node(const int node);
void setup_layout();
int node_number(int x, int y, int z, int t);
int node_index(int x, int y, int z, int t);
size_t num_sites(int node);
const int *get_logical_dimensions();
const int *get_logical_coordinate();
void get_coords(int coords[], int node, int index);

/* make_lattice.c */
void make_lattice();
void free_lattice();

/* nersc_cksum.c */
u_int32type nersc_cksum();

/* make_global_fields.c */
void make_global_fields();

/* plaquette.c */
void plaquette(Real *ss_plaq, Real *st_plaq);

// ploop.c
complex ploop(int dir);

// print_var3.c
void print_var3();

/* ranmom.c */
void ranmom();

/* remap standard I/O */
int remap_stdio_from_args(int argc, char *argv[]);

/* ranstuff.c */
void initialize_prn(double_prn *prn_pt, int seed, int index);
Real myrand(double_prn *prn_pt);

/* reunitarize.c */
int check_deviation();
void reunitarize();
// Use LAPACK singular value decomposition for reunitarization
// http://www.netlib.org/lapack/explore-3.1.1-html/zgesvd.f.html
// First and second arguments tell LAPACK to compute all singular values
// Third and fourth arguments are the dimensions of the matrix (both NCOL)
// Fifth argument is the input matrix (lost)
// Sixth argument is the leading dimension (NCOL)
// Seventh argument is the array of singular values (discarded)
// Eight argument is the matrix of left singular vectors (left)
// Ninth argument is the dimension of left (NCOL)
// Tenth argument is the matrix of right singular vectors (right^dag)
// Eleventh argument is the dimension of right^dag (NCOL)
// Twelfth argument is complex workspace of size given by the 13th argument
// Fourteenth argument is real workspace of size 5 * NCOL
// Final argument reports success or information about failure
void zgesvd_(char *A1, char *A2, int *N1, int *N2, double *store,
             int *lda, double *junk, double *left, int *Nl,
             double *right, int *Nr, double *work, int *Nwork,
             double *Rwork, int *stat);
#endif
// -----------------------------------------------------------------
