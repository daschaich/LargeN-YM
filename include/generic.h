// -----------------------------------------------------------------
// Macros and declarations for miscellaneous generic routines
#ifndef _GENERIC_H
#define _GENERIC_H

// Other generic directory declarations are elsewhere:
//   See comdefs.h for communications
//   See io_lat.h for I/O
//   See io_wprop.h for propagators
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
/* ape_smear.c */
void ape_smear_field(
  matrix *src,       /* Gauge field input unsmeared */
  matrix *dest,      /* Gauge field output smeared */
  Real staple_weight,    /* single staple weight */
  Real link_u0,          /* single link weight - used in normalization
                             if SU(3) projection is turned off */
  int space_only,         /* = 1 (true) smear space-like links with
                only spacelike staples
           = 0 (false) smear all links with
           all staples */
  int nhits,              /* reproject onto SU(3): number of
           SU(2) hits. 0 for no reprojection */
  Real tol               /* tolerance for SU(3) projection.
           If nonzero, treat nhits as a maximum
           number of hits.  If zero, treat nhits
           as a prescribed number of hits. */
  );

void ape_smear_dir(
  field_offset src,       /* field offset for matrix[4] type
           input unsmeared links */
  int dir1,               /* link direction to smear */
  field_offset dest,      /* field offset for matrix type
           pointing to a specific direction
           output smeared links */
  Real staple_weight,    /* single staple weight */
  Real link_u0,          /* single link weight - used in normalization
                             if SU(3) projection is turned off */
  int space_only,         /* = 1 (true) smear space-like links with
                only spacelike staples
           = 0 (false) smear all links with
           all staples */
  int nhits,              /* reproject onto SU(3): number of
           SU(2) hits. 0 for no reprojection */
  Real tol               /* tolerance for SU(3) projection.
           If nonzero, treat nhits as a maximum
           number of hits.  If zero, treat nhits
           as a prescribed number of hits. */
  );

void ape_smear(
  field_offset src,       /* field offset for matrix type
           input unsmeared links */
  field_offset dest,      /* field offset for matrix type
           output smeared links */
  Real staple_weight,    /* single staple weight */
  Real link_u0,          /* single link weight - used in normalization
                             if SU(3) projection is turned off */
  int space_only,         /* = 1 (true) smear space-like links with
                only spacelike staples
           = 0 (false) smear all links with
           all staples */
  int nhits,              /* reproject onto SU(3): number of
           SU(2) hits. 0 for no reprojection */
  Real tol               /* tolerance for SU(3) projection.
           If nonzero, treat nhits as a maximum
           number of hits.  If zero, treat nhits
           as a prescribed number of hits. */
  );

/* ax_gauge.c */
void ax_gauge();

/* check_unitarity.c */
Real check_unitarity();

// nersc_cksum.c
void linktrsum(double_complex *linktr);

// plaquette.c
void plaquette(double *ss_plaq, double *st_plaq);
void plaquette_frep(double *ss_plaq_frep, double *st_plaq_frep);

// plaquette_lcl.c
void plaquette_lcl(double *ss_plaq, double *st_plaq);
void plaquette_frep_lcl(double *ss_plaq_frep, double *st_plaq_frep);

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

/* reunitarize2.c */
void reunitarize();
#endif
// -----------------------------------------------------------------
