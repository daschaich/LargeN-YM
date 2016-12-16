#ifndef _GENERIC_H
#define _GENERIC_H
/************************ generic.h *************************************
*									*
*  Macros and declarations for miscellaneous generic routines           *
*  This header is for codes that call generic routines                  *
*  MIMD version 7 							*
*									*
*/

/* Other generic directory declarations are elsewhere:

   For com_*.c, see comdefs.h
   For io_lat4.c io_ansi.c, io_nonansi.c, io_piofs.c, io_romio.c see io_lat.h
   For io_prop_w.c, see io_wprop.h
*/

#include <stdio.h>
#include "../include/int32type.h"
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/macros.h"
#include "../include/random.h"
#include "../include/file_types.h"
#include "../include/io_lat.h"

/* ape_smear.c */

void ape_smear_field(
  su3_matrix *src,       /* Gauge field input unsmeared */
  su3_matrix *dest,      /* Gauge field output smeared */
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
  field_offset src,       /* field offset for su3_matrix[4] type
			     input unsmeared links */
  int dir1,               /* link direction to smear */
  field_offset dest,      /* field offset for su3_matrix type
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
  field_offset src,       /* field offset for su3_matrix type
			     input unsmeared links */
  field_offset dest,      /* field offset for su3_matrix type
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

/* bsd_sum.c */
int32type bsd_sum (char *data,int32type total_bytes);

/* check_unitarity.c */
Real check_unitarity( void );

/* d_linktrsum */
void d_linktrsum(double_complex *linktrsum);

/* d_plaq?.c */
void d_plaquette(double *ss_plaq, double *st_plaq);
void d_plaquette_frep(double *ss_plaq_frep, double *st_plaq_frep);
void d_plaquette_lcl(double *ss_plaq, double *st_plaq);
void d_plaquette_frep_lcl(double *ss_plaq_frep, double *st_plaq_frep);

/* field_strength.c */
void make_field_strength(
  field_offset link_src,       /* field offset for su3_matrix[4] type
				  for the source link matrices */
  field_offset field_dest      /* field offset for su3_matrix[6] type
				  for the resulting field strength */
  );

/* gaugefix.c and gaugefix2.c */
void gaugefix(int gauge_dir,Real relax_boost,int max_gauge_iter,
	      Real gauge_fix_tol, field_offset diffmat, field_offset sumvec,
	      int nvector, field_offset vector_offset[], int vector_parity[],
	      int nantiherm, field_offset antiherm_offset[],
	      int antiherm_parity[] );

/* gauge_force_imp.c and gauge_force_symzk1_qop.c */
void imp_gauge_force( Real eps, field_offset mom_off );

/* gauge_stuff.c */
double imp_gauge_action();
void g_measure();
void make_loop_table();
void dsdu_qhb_subl(int dir, int subl);
int get_max_length();
int get_nloop();
int get_nreps();
int *get_loop_length();
int *get_loop_num();
int ***get_loop_table();
Real **get_loop_coeff();

/* glueball_op.c */
void make_glueball_ops();
void measure_glueball_ops();

/* hvy_pot.c */
void hvy_pot( field_offset links );

/* io_detect.c */
int io_detect(char *filename, file_type ft[], int ntypes);

/* io_helpers.c */
gauge_file *save_lattice( int flag, char *filename, char *stringLFN );
gauge_file *reload_lattice( int flag, char *filename);
int ask_starting_lattice( FILE *fp, int prompt, int *flag, char *filename );
int ask_ending_lattice( FILE *fp, int prompt, int *flag, char *filename );
int ask_ildg_LFN(FILE *fp, int prompt, int flag, char *stringLFN);
void coldlat();
void funnylat();
#ifdef SF
int ask_sf_flag( FILE *fp, int prompt, char *flag_name_string, int *value );
#endif
int get_f( FILE *fp, int prompt, char *variable_name_string, Real *value );
int get_i( FILE *fp, int prompt, char *variable_name_string, int *value );
int get_vi( FILE *fp, int prompt, char *variable_name_string,
	    int *value, int nvalues );
int get_vf( FILE *fp, int prompt, char *variable_name_string,
	    Real *value, int nvalues );
int get_s( FILE *fp, int prompt, char *variable_name_string, char *value );
int get_prompt( FILE *fp, int *value );

/* layout_*.c */
void setup_layout( void );
int node_number(int x,int y,int z,int t);
int node_index(int x,int y,int z,int t);
size_t num_sites(int node);
const int *get_logical_dimensions();
const int *get_logical_coordinate();

/* make_lattice.c */
void make_lattice();
void free_lattice();

/* nersc_cksum.c */
u_int32type nersc_cksum( void );

/* make_global_fields.c */
void make_global_fields();

/* path_product.c */
void path_product( const int *dir, const int length, su3_matrix *tempmat1);
void path_prod_subl(const int *dir, const int length, const int subl,
		    su3_matrix *tempmat1);

/* plaquette4.c */
void plaquette(Real *ss_plaq,Real *st_plaq);

/* ploop?.c */
complex ploop( void );

/* ploop_staple.c */
complex ploop_staple(Real alpha_fuzz);

/* project_su3_hit.c */
void project_su3(
   su3_matrix *w,         /* input initial guess. output resulting
                             SU(3) matrix */
   su3_matrix *q,         /* starting 3 x 3 complex matrix */
   int Nhit,              /* number of SU(2) hits. 0 for no projection */
   Real tol              /* tolerance for SU(3) projection.
			     If nonzero, treat Nhit as a maximum
			     number of hits.  If zero, treat Nhit
			     as a prescribed number of hits. */
   );

/* rand_gauge.c */
void rand_gauge(field_offset G);

/* ranmom.c */
void ranmom( void );

/* remap standard I/O */
int remap_stdio_from_args(int argc, char *argv[]);

/* ranstuff.c */
void initialize_prn(double_prn *prn_pt, int seed, int index);
Real myrand(double_prn *prn_pt);

/* restrict_fourier.c */
void setup_restrict_fourier( int *key, int *slice);
void restrict_fourier(
     field_offset src,	 /* src is field to be transformed */
     field_offset space, /* space is working space, same size as src */
     field_offset space2,/* space2 is working space, same size as src */
                         /* space2 is needed only for non power of 2 */
     int size,		 /* Size of field in bytes.  The field must
			    consist of size/sizeof(complex) consecutive
			    complex numbers.  For example, an su3_vector
			    is 3 complex numbers. */
     int isign);	 /* 1 for x -> k, -1 for k -> x */

/* reunitarize2.c */
void reunitarize( void );
/* int reunit_su3(su3_matrix *c); [now declared only in source file generic/reunitarize2.c] */

/* show_generic_opts.c */
void show_generic_opts( void );

#ifdef QCDOC
void *qcdoc_alloc(size_t nbytes);
void qfree(void *);
#endif

#endif	/* _GENERIC_H */
