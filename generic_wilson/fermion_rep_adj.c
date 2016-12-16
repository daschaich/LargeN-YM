/*
 Create link field 'link' in rep of fermions from field 'linkf'
 in adjoint rep
*/
#include "generic_wilson_includes.h"

#if (DIMF != NCOL*NCOL-1)
	#error "Wrong version of fermion_rep!"
#endif
#if (FREP != adjoint)
        #error "Wrong version of fermion_rep!"
#endif


void make_fermion_rep_matrix(su3_matrix_f *a, su3_matrix *b){
int i,j,k,l,m,n,count;
float norm;

su3_matrix_f lambda[DIMF];
su3_matrix_f t1,t2;
su3_matrix_f square;
complex trace;


/*
	dumpmat_f(a);
*/


   for(i=0;i<DIMF;i++){
	   for(m=0;m<NCOL;m++)for(n=0;n<NCOL;n++)
	lambda[i].e[m][n]=cmplx(0.0,0.0);
   }

        count=0;
	for(l=2;l<=NCOL;l++){
	for(k=1;k<l;k++){
	lambda[count].e[k-1][l-1].real=lambda[count].e[l-1][k-1].real=0.5;
	lambda[count+1].e[k-1][l-1].imag=0.5;
	lambda[count+1].e[l-1][k-1].imag= -0.5;
	count+=2;
	}
	}

	for(i=0;i<NCOL-1;i++){
	k=i+1;
	norm=sqrt((float)(2*k*(k+1)));
	  for(j=1;j<=k;j++) lambda[i+count].e[j-1][j-1].real=1.0/norm;
	  lambda[i+count].e[k][k].real = (-(float)k)/norm;
	}

        for(i=0;i<DIMF;i++){
        printf("\n i= %d\n",i);
        dumpmat_f(&lambda[i]);
        mult_su3_nn_f(&lambda[i],&lambda[i],&square);
        trace=trace_su3_f(&square);
        printf("Trace tsq %e %e\n",trace.real,trace.imag);
        }

/*
	for(i=0;i<DIMF;i++)for(j=0;j<=i;j++){
	printf("\n i= %d  j=%d\n",i,j);fflush(stdout);
	mult_su3_nn_f(&lambda[i],&lambda[j],&square);
	trace=trace_su3_f(&square);
	printf("Trace tsq %e %e\n",trace.real,trace.imag);
	}
*/
	/* now b[i][j]= 0.5*Tr \lambda_i a^\dagger  \lambda[j] a */

	for(i=0;i<DIMF;i++)for(j=0;j<DIMF;j++){
	  	mult_su3_nn_f(&(lambda[j]),a,&t1);
	  	mult_su3_an_f(a,&t1,&t2);
		mult_su3_nn_f(&lambda[i],&t2,&t1);
		trace=trace_su3_f(&t1);
                CMULREAL(trace,2.0,b->e[i][j]);
	}

/*
	dumpmat(b);
*/

}	/* END make_fermion_rep_matrix() */




/*** make_fermion_rep_driv *******************************************
YS Sept 07

Generates parametric derivatives of boundary links using
the Leibniz rule, for use in measuring 1/g^2(L).

Minimilly adapted from make_fermion_rep_matrix.
The code assumes that the boundary matrices are diagonal

*/

#ifdef SF

void make_fermion_rep_driv(su3_matrix_f *a, su3_matrix *b,
                           Real *factor, int sign) {
  int i,j;
  /*

  for(i=0;i<NCOL;i++)for(j=0;j<NCOL;j++){
    b->e[i][j]=cmplx(0.0,0.0);
  }

  for(i=0;i<NCOL;i++) {
    b->e[i][i].imag = (Real)sign * factor[i] * a->e[i][i].real;
    b->e[i][i].real = -(Real)sign * factor[i] * a->e[i][i].imag;
  }
  */
  node0_printf("can't do make_fermion_rep_driv for adjoint...\n"); exit(1);
}	/* END make_fermion_rep_driv */

#endif /* SF */
