/*
Compute SUM_dirs (
    (1 + isign*gamma[dir]) * U(x, dir) * src(x+dir)
  + (1 - isign*gamma[dir]) * U_adj(x-dir, dir) * src(x-dir)
)
*/

#include "generic_wilson_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void dslash_w_site(field_offset src, field_offset dest, int isign,
                   int parity) {

  register int i, dir, otherparity = 0;
  register site *s;
  half_wilson_vector hwvx, hwvy, hwvz, hwvt;
  msg_tag *tag[8];

  switch(parity) {
    case EVEN:       otherparity = ODD;        break;
    case ODD:        otherparity = EVEN;       break;
    case EVENANDODD: otherparity = EVENANDODD; break;
  }

#ifdef MAXHTMP
  /* NOTE: We should be defining MAXHTMP in all applications using
     dslash and dslash_w_site */
  if (MAXHTMP < 8) {
    printf("dslash: MAXHTMP must be 8 or more!\n");
    terminate(1);
  }
#endif
  if (N_POINTERS < 8) {
    printf("dslash: N_POINTERS must be 8 or more!\n");
    terminate(1);
  }

  /* Take Wilson projection for src displaced in up direction, gather
     it to "our site" */
  FORSOMEPARITY(i, s, otherparity) {
    wp_shrink_4dir((wilson_vector *)F_PT(s, src), &(s->htmp[XUP]),
        &(s->htmp[YUP]), &(s->htmp[ZUP]), &(s->htmp[TUP]), isign);
  }
  FORALLUPDIR(dir) {
    tag[dir] = start_gather_site(F_OFFSET(htmp[dir]), sizeof(half_wilson_vector),
        dir, parity, gen_pt[dir]);
  }

  /* Take Wilson projection for src displaced in down direction,
     multiply it by adjoint link matrix, gather it "up" */
  FORSOMEPARITY(i, s, otherparity) {
    wp_shrink_4dir((wilson_vector *)F_PT(s, src),
        &hwvx, &hwvy, &hwvz, &hwvt, -isign);
    mult_adj_su3_mat_hwvec(&(s->link[XUP]), &hwvx, &(s->htmp[XDOWN]));
    mult_adj_su3_mat_hwvec(&(s->link[YUP]), &hwvy, &(s->htmp[YDOWN]));
    mult_adj_su3_mat_hwvec(&(s->link[ZUP]), &hwvz, &(s->htmp[ZDOWN]));
    mult_adj_su3_mat_hwvec(&(s->link[TUP]), &hwvt, &(s->htmp[TDOWN]));
  }

  FORALLUPDIR(dir) {
    tag[OPP_DIR(dir)] = start_gather_site(F_OFFSET(htmp[OPP_DIR(dir)]),
        sizeof(half_wilson_vector), OPP_DIR(dir),
        parity, gen_pt[OPP_DIR(dir)]);
  }

  /* Set dest to zero */
  /* Take Wilson projection for src displaced in up direction, gathered,
     multiply it by link matrix, expand it, and add.
     to dest */
  FORALLUPDIR(dir)
    wait_gather(tag[dir]);

  FORSOMEPARITY(i, s, parity) {
    mult_su3_mat_hwvec(&(s->link[XUP]),
        (half_wilson_vector *)(gen_pt[XUP][i]), &hwvx);
    mult_su3_mat_hwvec(&(s->link[YUP]),
        (half_wilson_vector *)(gen_pt[YUP][i]), &hwvy);
    mult_su3_mat_hwvec(&(s->link[ZUP]),
        (half_wilson_vector *)(gen_pt[ZUP][i]), &hwvz);
    mult_su3_mat_hwvec(&(s->link[TUP]),
        (half_wilson_vector *)(gen_pt[TUP][i]), &hwvt);
    grow_add_four_wvecs((wilson_vector *)F_PT(s, dest),
        &hwvx, &hwvy, &hwvz, &hwvt, isign);
  }
  FORALLUPDIR(dir)
    cleanup_gather(tag[dir]);

  /* Take Wilson projection for src displaced in down direction,
     expand it, and add to dest */
  FORALLUPDIR(dir)
    wait_gather(tag[OPP_DIR(dir)]);

  FORSOMEPARITY(i, s, parity) {
    grow_sum_four_wvecs((wilson_vector *)F_PT(s, dest),
        (half_wilson_vector *)(gen_pt[XDOWN][i]),
        (half_wilson_vector *)(gen_pt[YDOWN][i]),
        (half_wilson_vector *)(gen_pt[ZDOWN][i]),
        (half_wilson_vector *)(gen_pt[TDOWN][i]),
        -isign);
  }
  FORALLUPDIR(dir)
    cleanup_gather(tag[OPP_DIR(dir)]);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
/* Special dslash for use by congrad.  Uses restart_gather_site() when
  possible. Last argument is an integer, which will tell if
  gathers have been started.  If is_started = 0, use
  start_gather_site, otherwise use restart_gather_site.
  Argument "tag" is a vector of a msg_tag *'s to use for
  the gathers.
  The calling program must clean up the gathers! */
void dslash_w_site_special(field_offset src, field_offset dest,
                           int isign, int parity, msg_tag **tag,
                           int is_started) {

  half_wilson_vector hwvx, hwvy, hwvz, hwvt;

  register int i;
  register site *s;
  register int dir, otherparity = 0;

  switch(parity) {
    case EVEN:       otherparity = ODD;        break;
    case ODD:        otherparity = EVEN;       break;
    case EVENANDODD: otherparity = EVENANDODD; break;
  }

#ifdef MAXHTMP
  /* NOTE: We should be defining MAXHTMP in all applications using
     dslash and dslash_w_site */
  if (MAXHTMP < 8) {
    printf("dslash_w_site_special: MAXHTMP must be 8 or more!\n");
    terminate(1);
  }
#endif
  if (N_POINTERS < 8) {
    printf("dslash_w_site_special: N_POINTERS must be 8 or more!\n");
    terminate(1);
  }

  /* Take Wilson projection for src displaced in up direction, gather
     it to "our site" */
  FORSOMEPARITY(i, s, otherparity) {
    wp_shrink_4dir((wilson_vector *)F_PT(s, src), &(s->htmp[XUP]),
        &(s->htmp[YUP]), &(s->htmp[ZUP]), &(s->htmp[TUP]), isign);
  }

  FORALLUPDIR(dir) {
    if (is_started == 0)
      tag[dir] = start_gather_site(F_OFFSET(htmp[dir]),
        sizeof(half_wilson_vector), dir, parity, gen_pt[dir]);
    else restart_gather_site(F_OFFSET(htmp[dir]),
        sizeof(half_wilson_vector), dir, parity, gen_pt[dir],
        tag[dir]);
  }

  /* Take Wilson projection for src displaced in down direction,
     multiply it by adjoint link matrix, gather it "up" */
  FORSOMEPARITY(i, s, otherparity) {
    wp_shrink_4dir((wilson_vector *)F_PT(s, src),
        &hwvx, &hwvy, &hwvz, &hwvt, -isign);
    mult_adj_su3_mat_hwvec(&(s->link[XUP]), &hwvx, &(s->htmp[XDOWN]));
    mult_adj_su3_mat_hwvec(&(s->link[YUP]), &hwvy, &(s->htmp[YDOWN]));
    mult_adj_su3_mat_hwvec(&(s->link[ZUP]), &hwvz, &(s->htmp[ZDOWN]));
    mult_adj_su3_mat_hwvec(&(s->link[TUP]), &hwvt, &(s->htmp[TDOWN]));
  }

  FORALLUPDIR(dir) {
    if (is_started == 0)
      tag[OPP_DIR(dir)] = start_gather_site(
        F_OFFSET(htmp[OPP_DIR(dir)]), sizeof(half_wilson_vector),
        OPP_DIR(dir), parity, gen_pt[OPP_DIR(dir)]);
    else
      restart_gather_site(
        F_OFFSET(htmp[OPP_DIR(dir)]), sizeof(half_wilson_vector),
        OPP_DIR(dir), parity, gen_pt[OPP_DIR(dir)],
        tag[OPP_DIR(dir)]);
  }

  /* Set dest to zero */
  /* Take Wilson projection for src displaced in up direction, gathered,
     multiply it by link matrix, expand it, and add.
     to dest */
  FORALLUPDIR(dir)
    wait_gather(tag[dir]);

  FORSOMEPARITY(i, s, parity) {
    mult_su3_mat_hwvec(&(s->link[XUP]),
        (half_wilson_vector *)(gen_pt[XUP][i]), &hwvx);
    mult_su3_mat_hwvec(&(s->link[YUP]),
        (half_wilson_vector *)(gen_pt[YUP][i]), &hwvy);
    mult_su3_mat_hwvec(&(s->link[ZUP]),
        (half_wilson_vector *)(gen_pt[ZUP][i]), &hwvz);
    mult_su3_mat_hwvec(&(s->link[TUP]),
        (half_wilson_vector *)(gen_pt[TUP][i]), &hwvt);
    grow_add_four_wvecs((wilson_vector *)F_PT(s, dest),
        &hwvx, &hwvy, &hwvz, &hwvt, isign);
  }

  /* Take Wilson projection for src displaced in down direction,
     expand it, and add to dest */
  FORALLUPDIR(dir)
    wait_gather(tag[OPP_DIR(dir)]);

  FORSOMEPARITY(i, s, parity) {
    grow_sum_four_wvecs((wilson_vector *)F_PT(s, dest),
        (half_wilson_vector *)(gen_pt[XDOWN][i]),
        (half_wilson_vector *)(gen_pt[YDOWN][i]),
        (half_wilson_vector *)(gen_pt[ZDOWN][i]),
        (half_wilson_vector *)(gen_pt[TDOWN][i]),
        -isign);
  }
}
// -----------------------------------------------------------------
