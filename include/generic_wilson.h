// -----------------------------------------------------------------
// Macros and declarations for generic Wilson routines
#ifndef _GENERIC_WILSON_H
#define _GENERIC_WILSON_H

#include "../include/su3.h"
#include "../include/comdefs.h"
#include "../include/macros.h"
#include "../include/generic_quark_types.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
/* For various inversion routines.  Not used sytematically yet */
enum guess_params { START_ZERO_GUESS = 0 ,  START_NONZERO_GUESS };

int wilson_invert(     // Return value is number of iterations taken
    field_offset src,   // type wilson_vector (source already created)
    field_offset dest,  // type wilson_vector (answer and initial guess)
    field_offset sav,   // type wilson_vector (for saving source)
    int (*invert_func)(field_offset src, field_offset dest,
      quark_invert_control *qic,
      void *dmp),
    quark_invert_control *qic, /* inverter control */
    void *dmp            /* Passthrough Dirac matrix parameters */
    );

void boundary_flip    (int sign);
void boundary_flip_x  (int sign);
void boundary_flip_z3 (int sign);

void copy_site_wilson_vector(field_offset src, field_offset dest);

/* For quark source routines */
/* The Weyl representation types are included for w_source_h */
enum source_type {
  POINT = 1, GAUSSIAN, CUTOFF_GAUSSIAN,
  POINT_WEYL, CUTOFF_GAUSSIAN_WEYL } ;

/* w_source.c */
void w_source(field_offset src,wilson_quark_source *wqs);
void w_sink(field_offset snk,wilson_quark_source *wqs);
void w_sink_scalar(field_offset snk,wilson_quark_source *wqs);
int ask_quark_source(int prompt, int *type, char *descrp );

/* w_source_h.c (also has an ask_quark_source) */
Real *make_template(Real gamma, int cutoff);
void w_source_h(field_offset src,wilson_quark_source *wqs);
void free_source_template();

void bj_to_weyl(wilson_vector *src, wilson_vector *dest);
void dslash_w_site(field_offset src,field_offset dest,
                   int isign, int parity);
void dslash_w_site_special(field_offset src,field_offset dest,
                           int isign, int parity,msg_tag **tag,
                           int is_started);
void dslash_w_field(wilson_vector *src, wilson_vector *dest,
                    int isign, int parity);
void dslash_w_field_special(wilson_vector *src, wilson_vector *dest,
                            int isign, int parity,msg_tag **tag,
                            int is_started);
void dslash_w_3D(field_offset src, field_offset dest, int isign, int parity);

void cleanup_dslash_temps();
void cleanup_tmp_links();
void meson_cont(field_offset src1,field_offset src2,
                int *gamma_in, int *gamma_out, int n_in, int n_out,
                complex *prop);

void w_meson(field_offset src1,field_offset src2,complex *prop[10]);
void d_w_meson(field_offset src1,field_offset src2,double_complex *prop[10]);
void w_baryon(field_offset src1,field_offset src2,field_offset src3,
              complex *prop[4]);

void w_baryon_hl(field_offset src1,field_offset src2,
     field_offset src3, complex *prop[6]);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Higher-rep routines
void fermion_rep();
void make_fermion_rep_matrix(su3_matrix_f *a, su3_matrix *b);
#endif
// -----------------------------------------------------------------
