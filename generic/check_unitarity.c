/*********************** check_unitarity.c ***************************/
/* MIMD version 7 */
/* Claude Bernard, original version */
/* Modified 7/96 DT to quit after finding violation */
/* 9/4/96 DT added Urs' row orthogonality checks, if STRONG defined */
/* 11/2/98 CD fixed return value of max deviation */
/* 01/20/00 UMH combined with Schroedinger functional version */

#define STRONG	/* check row orthogonality as well as norms */

/* Check unitarity of link matrices, quit if not unitary */

#include "generic_includes.h"

#define TOLERANCE (0.0001)
/*#define UNIDEBUG */
Real check_su3(su3_matrix_f *c);

Real check_unitarity() {
register int i,dir;
int ii,jj;
register site *s;
register su3_matrix_f *mat;
Real deviation,max_deviation;
double av_deviation;
 union {
   Real fval;
   int ival;
 } ifval;
 
 max_deviation=av_deviation=0;
 FORALLSITES(i,s){
#ifdef SCHROED_FUN
   for(dir=XUP; dir<=TUP; dir++ ) if(dir==TUP || s->t>0){
#else
   for(dir=XUP; dir<=TUP; dir++ ){
#endif
     mat = (su3_matrix_f *)&(s->linkf[dir]);
     deviation=check_su3( mat );
     if (deviation>TOLERANCE){
       printf("Unitarity problem on node %d, site %d, dir %d, deviation=%f\n",
	      mynode(),i,dir,deviation);
       printf("SU(N) matrix:\n");
       for(ii=0;ii<NCOL;ii++){
	 for(jj=0;jj<NCOL;jj++){
	   printf("%f ",(*mat).e[ii][jj].real); 
	   printf("%f ",(*mat).e[ii][jj].imag); 
	 }
	 printf("\n");
       }
       printf("repeat in hex:\n");
       for(ii=0;ii<NCOL;ii++){
	 for(jj=0;jj<NCOL;jj++){
	   ifval.fval = (*mat).e[ii][jj].real; 
	   printf("%08x ", ifval.ival); 
	   ifval.fval = (*mat).e[ii][jj].imag; 
	   printf("%08x ", ifval.ival); 
	 }
	 printf("\n");
       }
       printf("  \n \n");
       fflush(stdout); terminate(1);
     }
     if(max_deviation<deviation) max_deviation=deviation;
     av_deviation += deviation*deviation;
   }
 }
 av_deviation = sqrt(av_deviation/(4*i));
#ifdef UNIDEBUG
 printf("Deviation from unitarity on node %d: max %g, avrg %g\n",
	mynode(), max_deviation, av_deviation);
#endif
 if(max_deviation> TOLERANCE) 
   printf("Unitarity problem on node %d, maximum deviation=%f\n",
	  mynode(),max_deviation);
 return max_deviation;
}  /*check_unitarity() */

Real check_su3(su3_matrix_f *c) {
  register Real ar, ai, ari, max;
  register int i,j,k;
  
  /* first check  normalization of each row */
  for(i=0 , max=0.; i<NCOL; ++i) {
    ar=0.0;
    for(j=0;j<NCOL;j++){
      ar += (*c).e[i][j].real * (*c).e[i][j].real +    /* sum of squares of row */
	(*c).e[i][j].imag * (*c).e[i][j].imag;
    }
    ar =  fabs( sqrt((double)ar) - 1.);
    if(max<ar) max=ar;
  }
#ifdef STRONG
  
  /* Test orthogonality of row i and row j */
  for(i=0;i<NCOL; ++i) {
    for(j=i+1;j<NCOL;j++){
      ar=0.0; ai=0.0;
      for(k=0;k<NCOL;k++){
	ar += (*c).e[i][k].real * (*c).e[j][k].real +     /* real part of 0 dot 1 */
	  (*c).e[i][k].imag * (*c).e[j][k].imag;
	ai += (*c).e[i][k].real * (*c).e[j][k].imag -     /* imag part of 0 dot 1 */
	  (*c).e[i][k].imag * (*c).e[j][k].real ;
      }
      ari = sqrt((double)(ar*ar + ai*ai));
      if(max<ari) max=ari;
    }
  }  
#endif /*STRONG*/
    
  return(max);
  
} /* check_su3 */
