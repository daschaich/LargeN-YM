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

void boundary_flip(int sign);

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

void cleanup_dslash_temps();
void cleanup_tmp_links();
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Higher-rep routines
void fermion_rep();
void make_fermion_rep_matrix(su3_matrix_f *a, su3_matrix *b);
#endif
// -----------------------------------------------------------------
