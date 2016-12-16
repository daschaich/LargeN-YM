/************************** w_source_h.c *****************************/
/* MIMD version 7 */
/* Set up a quark source based on a Gaussian with a box cutoff */

/* CB wrote this */
/* 11/25/97 modifications for version 5 CD */
/*  2/14/98 corrected memory leak CD */
/*  2/15/98 to fail on illegal source type CD */

/* Initialize a source for the inverter */
#include "generic_wilson_includes.h"
#include <string.h>

static Real *source_wall_template = NULL;

/* Make precomputed Gaussian weights */

Real *make_template(Real gamma, int cutoff)  
{
  
  int rx,ry,rz;
  Real radius2,scale;
  Real *my_template;

  /* Gaussian trial source centered on  x0,y0,z0 in each timeslice t0 */
  
  my_template =(Real *)malloc( 
				 (cutoff+1)*(cutoff+1)*(cutoff+1)*sizeof(Real) );
  if(my_template==NULL){
    printf("NODE %d: no room for source_wall_template\n",this_node);
    terminate(1);
  }
  
  for(rz=0;rz<=cutoff;rz++)for(ry=0;ry<=cutoff;ry++)for(rx=0;rx<=cutoff;rx++){
    radius2= (Real)(rx*rx + ry*ry + rz*rz);
    scale = (Real)exp((double)(- radius2*gamma));  
    my_template[rx + (cutoff+1)*(ry +(cutoff+1)*rz)]=scale;
  }
  
  return my_template;
} /* make_template */

void w_source_h(field_offset src,wilson_quark_source *wqs)
{
  register int i;
  register site *s; 
  Real size_src;
  
  int rx,ry,rz;
  Real scale;
  int x0,y0,z0,t0;
  int color,spin;
  int source_type,parity;
  int wall_cutoff;
  Real r0;
  
  wilson_vector  bj, weyl;
  
  /* Unpack structure */
  x0 = wqs->x0; y0 = wqs->y0; z0 = wqs->z0; t0 = wqs->t0;
  color = wqs->color; spin = wqs->spin;
  source_type = wqs->type;
  r0 = wqs->r0;
  parity= wqs->parity;
  wall_cutoff = wqs->wall_cutoff;
  
  /* Make template if not already done -- needed only for Gaussian sources */
  /* Make template if not already done */
  if(source_wall_template == NULL && 
     ((source_type == CUTOFF_GAUSSIAN) || 
      (source_type == CUTOFF_GAUSSIAN_WEYL)))
    source_wall_template = make_template(1./(r0*r0),wall_cutoff);
  
  /* zero src to be safe */
  FORALLSITES(i,s) {
    clear_wvec((wilson_vector *)F_PT(s,src)); 
  }
  /* Normalisation  */
  size_src=0.0;
  FORALLSITES(i,s) {
    size_src += magsq_wvec( ((wilson_vector *)F_PT(s,src)) );
  }
  g_floatsum(&size_src);
  size_src = (Real)sqrt((double)size_src);
  /*if(this_node==0)printf("WSOURCE 1--size_src=%e\n",
    (double)size_src); */
  
  /* make source vector in B&D convention & convert to Weyl conv. */
  clear_wvec( &bj );
  bj.d[spin].c[color].real = 1.0;
  
  if(source_type == POINT_WEYL || source_type == CUTOFF_GAUSSIAN_WEYL)
    {
      bj_to_weyl( &bj, &weyl); 
      node0_printf("w_source_h:: source rep. flipped bj-->weyl\n"); 
    }
  else
    {
      weyl = bj ; 
    }

  /*printf(" spin= %d, color= %d, weyl source: \n",spin,color);
    for(j=0;j<3;j++){
    printf(" weyl color= %d\n",j);
    for(k=0;k<4;k++){
    printf("\t (%e, %e)\n",weyl.d[k].c[j].real,
    weyl.d[k].c[j].imag);
    }
    }
    fflush(stdout);  */
  
  
  if(source_type == POINT || 
     source_type == POINT_WEYL) {
    /* load 1.0 into source at cooordinates given by source_coord */
    /* initialize src to be a delta function at point x0,y0,z0,t0 */
    
    if(node_number(x0,y0,z0,t0) == mynode()){
      i = node_index(x0,y0,z0,t0);
      *((wilson_vector *)F_PT(&(lattice[i]),src)) = weyl;
    }
  }
  else if(source_type == CUTOFF_GAUSSIAN || 
	  source_type == CUTOFF_GAUSSIAN_WEYL) {
    /* Formerly called "WALL" */
    /* Gaussian trial source of fixed parity 
       points centered on  x0,y0,z0,t0;  
       cutoff a rectangular distance wall_cutoff from the source */
    
    FORSOMEPARITY(i,s,parity) {
      if(s->t != t0) continue;	/* only do this if t==t0 */
      rz=abs(s->z - z0);
      if(rz > (nz/2) ) rz = nz - rz;
      if(rz>wall_cutoff) continue; /* don't do anything outside wall_cutoff */
      ry=abs(s->y - y0);
      if(ry > (ny/2) ) ry = ny - ry;
      if(ry>wall_cutoff) continue; /* don't do anything outside wall_cutoff */
      rx=abs(s->x - x0);
      if(rx > (nx/2) ) rx = nx - rx;
      if(rx>wall_cutoff) continue; /* don't do anything outside wall_cutoff */
      scale = source_wall_template[rx + (wall_cutoff+1)*(ry +(wall_cutoff+1)*rz)];
      scalar_mult_wvec( &weyl, scale, ((wilson_vector *)F_PT(s,src)) );  
    }
  }
  else {
    node0_printf("w_source_h: Unrecognized source type %d\n",source_type);
    terminate(1);
  }
  /* Normalisation  */
  size_src=0.0;
  FORALLSITES(i,s) {
    size_src += magsq_wvec( ((wilson_vector *)F_PT(s,src)) );
  }
  g_floatsum(&size_src);
  size_src = (Real)sqrt((double)size_src);
  /*if(this_node==0)printf("WSOURCE 2--size_src=%e\n",
    (double)size_src);  */
  
} /* w_source_h */

int ask_quark_source( int prompt, int *source_type, char *descrp)
{
  char savebuf[256];
  int status;
  
  if (prompt!=0)
    printf("enter 'point', or 'gaussian' for source type\n");
  status = scanf("%s",savebuf);
  if(status !=1) {
    printf("ask_quark_source: ERROR IN INPUT: source type command\n");
    return 1;
  }
  if(strcmp("point",savebuf) == 0 ){
    *source_type = POINT;
    strcpy(descrp,"point");
  }
  else if(strcmp("cutoff_gaussian",savebuf) == 0 ) {
    *source_type = CUTOFF_GAUSSIAN;
    descrp = "cutoff_gaussian";
  }
  else if(strcmp("point_weyl",savebuf) == 0 ) {
    *source_type = POINT_WEYL;
    strcpy(descrp,"point_weyl");
  }
  else if(strcmp("cutoff_gaussian_weyl",savebuf) == 0 ) {
    *source_type = CUTOFF_GAUSSIAN_WEYL;
    strcpy(descrp,"cutoff_gaussian_weyl");
  }
  else{
    printf("ask_source: ERROR IN INPUT: source command is invalid\n"); 
    return 1;
  }

  printf("%s\n",savebuf);
  return 0;
} /* ask_quark_source */

/* Free wall template malloc */
void free_source_template()
{
  if(source_wall_template != NULL)free(source_wall_template);
  source_wall_template = NULL;
}

/* So far no one has written the corresponding w_sink_h routine 
   11/25/97 CD */
