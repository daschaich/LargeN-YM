/*************  dslash_w_site.c *******************************/
/* MIMD version 7 */

/* 9/06/03 added gathers from temp - CD */

/*
	dslash_w_site(F_OFFSET(psi),F_OFFSET(mp),isign,l_parity);
Compute SUM_dirs (
    ( 1 + isign*gamma[dir] ) * U(x,dir) * src(x+dir)
  + ( 1 - isign*gamma[dir] ) * U_adj(x-dir,dir) * src(x-dir)
)

*/

/***  handling SF  *********************************************

In set_boundary_values() all wilson vector at t=0 are set to zero.
We gather from, but never to, t=0="nt" .
As a result, all wilson vectors at t=0 always remain zero,
and any gather from t=0 will always add zero on t=1 or t=nt-1.
*/

#include "generic_wilson_includes.h"
/* Temporary work space for dslash_w_field and dslash_w_field_special */
static half_wilson_vector *htmp[8] ;

/* Flag indicating if temp is allocated               */
static int temp_not_allocated=1 ;

/* Temporary work space for gauge links */
static su3_matrix *t_links;

/* Flag indicating if temp is allocated               */
static int tmp_links_not_set = 1;

void malloc_dslash_temps(){
  int j;

  if(!temp_not_allocated)return;
  for( j = 0; j < 8; j++ ){
    htmp[j] =(half_wilson_vector *)malloc(sites_on_node*sizeof(half_wilson_vector));
    if(htmp[j] == NULL){
      printf("node %d can't malloc htmp[%d]\n",this_node,j);
      terminate(1);
    }
  }
  temp_not_allocated = 0 ;
}

void cleanup_dslash_temps(){
  register int i ;
  if(!temp_not_allocated)
    for(i=0;i<8;i++) {
      free(htmp[i]) ;
    }
  temp_not_allocated=1 ;
}

void setup_tmp_links(){
  register int i, dir;
  register site *s;

  if(!tmp_links_not_set)return;
  t_links =(su3_matrix *)malloc(sites_on_node*4*sizeof(su3_matrix));
  if(t_links == NULL){
      printf("node %d can't malloc t_links\n",this_node);
      terminate(1);
  }
  FORALLSITES(i,s){
    FORALLUPDIR(dir){
      t_links[4*i+dir] = lattice[i].link[dir];
    }
  }
  tmp_links_not_set = 0 ;
}

void cleanup_tmp_links(){
  if(!tmp_links_not_set)
      free(t_links) ;
  tmp_links_not_set=1 ;
}


void dslash_w_field( wilson_vector *src, wilson_vector *dest, int isign, int parity) {
half_wilson_vector hwvx,hwvy,hwvz,hwvt;

register int i;
register site *s;
register int dir,otherparity=0;
msg_tag *tag[8];
su3_matrix *linkx,*linky,*linkz,*linkt;

  /* The calling program must clean up the temps! */
  malloc_dslash_temps();
  setup_tmp_links();

    switch(parity) {
	case EVEN:      otherparity=ODD; break;
	case ODD:       otherparity=EVEN; break;
	case EVENANDODD:        otherparity=EVENANDODD; break;
    }

    if(N_POINTERS < 8){
      printf("dslash: N_POINTERS must be 8 or more!\n");
      terminate(1);
     }


    /* Take Wilson projection for src displaced in up direction, gather
       it to "our site" */
    FORSOMEPARITY(i,s,otherparity){
        wp_shrink_4dir( &(src[i]), &(htmp[XUP][i]),
	    &(htmp[YUP][i]), &(htmp[ZUP][i]), &(htmp[TUP][i]), isign);
    }
    for( dir=XUP; dir <= TUP; dir++) {
	tag[dir]=start_gather_field( htmp[dir], sizeof(half_wilson_vector),
	    dir, parity, gen_pt[dir] );
    }

        /* Take Wilson projection for src displaced in down direction,
        multiply it by adjoint link matrix, gather it "up" */
    FORSOMEPARITY(i,s,otherparity){
        wp_shrink_4dir( &(src[i]), &hwvx, &hwvy, &hwvz, &hwvt, -isign);

	linkx = &t_links[4*i+XUP];
	linky = &t_links[4*i+YUP];
	linkz = &t_links[4*i+ZUP];
	linkt = &t_links[4*i+TUP];
	mult_adj_su3_mat_hwvec( linkx, &hwvx, &(htmp[XDOWN][i]));
	mult_adj_su3_mat_hwvec( linky, &hwvy, &(htmp[YDOWN][i]));
	mult_adj_su3_mat_hwvec( linkz, &hwvz, &(htmp[ZDOWN][i]));
	mult_adj_su3_mat_hwvec( linkt, &hwvt, &(htmp[TDOWN][i]));
    }

    for( dir=XUP; dir <= TUP; dir++) {
	tag[OPP_DIR(dir)]=start_gather_field(htmp[OPP_DIR(dir)],
		sizeof(half_wilson_vector), OPP_DIR(dir),
		parity, gen_pt[OPP_DIR(dir)] );
    }


	/* Set dest to zero */
        /* Take Wilson projection for src displaced in up direction, gathered,
		multiply it by link matrix, expand it, and add.
		to dest */
    for( dir=XUP; dir <= TUP; dir++) {
	wait_gather(tag[dir]);
    }
    FORSOMEPARITYDOMAIN(i,s,parity){
      linkx = &t_links[4*i+XUP];
      linky = &t_links[4*i+YUP];
      linkz = &t_links[4*i+ZUP];
      linkt = &t_links[4*i+TUP];
      mult_su3_mat_hwvec( linkx,
			   (half_wilson_vector * )(gen_pt[XUP][i]), &hwvx );
      mult_su3_mat_hwvec( linky,
			   (half_wilson_vector * )(gen_pt[YUP][i]), &hwvy );
      mult_su3_mat_hwvec( linkz,
			   (half_wilson_vector * )(gen_pt[ZUP][i]), &hwvz );
      mult_su3_mat_hwvec( linkt,
			   (half_wilson_vector * )(gen_pt[TUP][i]), &hwvt );
      grow_add_four_wvecs( &(dest[i]),
	    &hwvx, &hwvy, &hwvz, &hwvt, isign, 0 ); /* "0" is NOSUM */
    }
    for( dir=XUP; dir <= TUP; dir++) {
	cleanup_gather(tag[dir]);
    }

        /* Take Wilson projection for src displaced in down direction,
        expand it, and add to dest */
    for( dir=XUP; dir <= TUP; dir++) {
	wait_gather(tag[OPP_DIR(dir)]);
    }

    FORSOMEPARITYDOMAIN(i,s,parity){
	grow_add_four_wvecs( &(dest[i]),
	    (half_wilson_vector *)(gen_pt[XDOWN][i]),
	    (half_wilson_vector *)(gen_pt[YDOWN][i]),
	    (half_wilson_vector *)(gen_pt[ZDOWN][i]),
	    (half_wilson_vector *)(gen_pt[TDOWN][i]),
	    -isign, 1 );	/* "1" SUMs in current dest */
    }
    for( dir=XUP; dir <= TUP; dir++) {
	cleanup_gather(tag[OPP_DIR(dir)]);
    }

} /* end (of dslash_w_field() ) */


/* Special dslash for use by congrad.  Uses restart_gather_temp() when
  possible. Last argument is an integer, which will tell if
  gathers have been started.  If is_started=0,use
  start_gather_field, otherwise use restart_gather_field.
  Argument "tag" is a vector of a msg_tag *'s to use for
  the gathers.
  The calling program must clean up the gathers! */
void dslash_w_field_special( wilson_vector *src, wilson_vector *dest,
    int isign,int parity,msg_tag **tag,int is_started)

{
half_wilson_vector hwvx,hwvy,hwvz,hwvt;

register int i;
register site *s;
register int dir,otherparity=0;
su3_matrix *linkx,*linky,*linkz,*linkt;
  /* allocate temporary work space only if not already allocated */
  /* The calling program must clean up this space */
  malloc_dslash_temps();
  setup_tmp_links();

    switch(parity) {
	case EVEN:      otherparity=ODD; break;
	case ODD:       otherparity=EVEN; break;
	case EVENANDODD:        otherparity=EVENANDODD; break;
    }

    if(N_POINTERS < 8){
      printf("dslash_w_field_special: N_POINTERS must be 8 or more!\n");
      terminate(1);
     }


    /* Take Wilson projection for src displaced in up direction, gather
       it to "our site" */
    FORSOMEPARITY(i,s,otherparity){
        wp_shrink_4dir( &(src[i]), &(htmp[XUP][i]),
	    &(htmp[YUP][i]), &(htmp[ZUP][i]), &(htmp[TUP][i]), isign);
    }

    for( dir=XUP; dir <= TUP; dir++) {
	if(is_started==0)tag[dir]=start_gather_field( htmp[dir],
	    sizeof(half_wilson_vector), dir, parity, gen_pt[dir] );
	else restart_gather_field( htmp[dir],
	    sizeof(half_wilson_vector), dir, parity, gen_pt[dir],
			     tag[dir] );
    }

        /* Take Wilson projection for src displaced in down direction,
        multiply it by adjoint link matrix, gather it "up" */
    FORSOMEPARITY(i,s,otherparity){
      linkx = &t_links[4*i+XUP];
      linky = &t_links[4*i+YUP];
      linkz = &t_links[4*i+ZUP];
      linkt = &t_links[4*i+TUP];
      wp_shrink_4dir( &(src[i]), &hwvx, &hwvy, &hwvz, &hwvt, -isign);
      mult_adj_su3_mat_hwvec( linkx, &hwvx, &(htmp[XDOWN][i]));
      mult_adj_su3_mat_hwvec( linky, &hwvy, &(htmp[YDOWN][i]));
      mult_adj_su3_mat_hwvec( linkz, &hwvz, &(htmp[ZDOWN][i]));
      mult_adj_su3_mat_hwvec( linkt, &hwvt, &(htmp[TDOWN][i]));
    }

    for( dir=XUP; dir <= TUP; dir++) {
	if(is_started==0)tag[OPP_DIR(dir)]=start_gather_field(
	    htmp[OPP_DIR(dir)], sizeof(half_wilson_vector),
	    OPP_DIR(dir), parity, gen_pt[OPP_DIR(dir)] );
	else restart_gather_field(
	    htmp[OPP_DIR(dir)], sizeof(half_wilson_vector),
	    OPP_DIR(dir), parity, gen_pt[OPP_DIR(dir)],
	    tag[OPP_DIR(dir)] );
    }

	/* Set dest to zero */
        /* Take Wilson projection for src displaced in up direction, gathered,
		multiply it by link matrix, expand it, and add.
		to dest */
    for( dir=XUP; dir <= TUP; dir++) {
	wait_gather(tag[dir]);
    }
    FORSOMEPARITYDOMAIN(i,s,parity){
      linkx = &t_links[4*i+XUP];
      linky = &t_links[4*i+YUP];
      linkz = &t_links[4*i+ZUP];
      linkt = &t_links[4*i+TUP];
      mult_su3_mat_hwvec( linkx,
			   (half_wilson_vector * )(gen_pt[XUP][i]), &hwvx );
      mult_su3_mat_hwvec( linky,
			   (half_wilson_vector * )(gen_pt[YUP][i]), &hwvy );
      mult_su3_mat_hwvec( linkz,
			   (half_wilson_vector * )(gen_pt[ZUP][i]), &hwvz );
      mult_su3_mat_hwvec( linkt,
			   (half_wilson_vector * )(gen_pt[TUP][i]), &hwvt );
      grow_add_four_wvecs( &(dest[i]),
	    &hwvx, &hwvy, &hwvz, &hwvt, isign, 0 ); /* "0" is NOSUM */
    }

        /* Take Wilson projection for src displaced in down direction,
        expand it, and add to dest */
    for( dir=XUP; dir <= TUP; dir++) {
	wait_gather(tag[OPP_DIR(dir)]);
    }

    FORSOMEPARITYDOMAIN(i,s,parity){
	grow_add_four_wvecs( &(dest[i]),
	    (half_wilson_vector *)(gen_pt[XDOWN][i]),
	    (half_wilson_vector *)(gen_pt[YDOWN][i]),
	    (half_wilson_vector *)(gen_pt[ZDOWN][i]),
	    (half_wilson_vector *)(gen_pt[TDOWN][i]),
	    -isign, 1 );	/* "1" SUMs in current dest */
    }

} /* end (of dslash_w_field_special() ) */

