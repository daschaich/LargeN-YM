/******************* make_clov2.c **************************************/

/* MIMD version 7 */
/* version of 1/4/95 by UMH */
/* 1/21/00 combined with Schroedinger functional version - UMH */

/* Prepare the "clover" term */

#include "generic_clover_includes.h"
#include "../include/prefetch.h"
#define FETCH_UP 1
#define LOOPEND
#include "../include/loopend.h"

/******************* global_clov *********************************/
/* The global clover term.  Used when we are working with only one
   clover term at a time */

static clover *global_clov;

/******************* create_clov ***************************************/
/* Allocate space for a clover term */

clover *create_clov(void){

  clover *my_clov = (clover *)malloc(sizeof(clover));

  my_clov->clov      = (triangular *)malloc(sites_on_node*sizeof(triangular));
  my_clov->clov_diag = (diagonal *)malloc(sites_on_node*sizeof(diagonal));

  if(my_clov->clov == NULL || my_clov->clov_diag ==NULL){
    printf("create_clov(%d): malloc failed\n",this_node);
    free(my_clov);
    return NULL;
  }
  else
    return my_clov;
}


/******************* compute_clov ***************************************/
/* Computes the clover field */
/* Modified for DIMF != 3  -bqs 12.06	*/

void compute_clov(clover *my_clov, Real Clov_c)
{

  triangular *clov = my_clov->clov;
  diagonal *clov_diag = my_clov->clov_diag;
  register int i;
  register site *s;

  register int j,k,jk,jk2;
  register complex ctmp;
  su3_matrix *f_mn;
  char myname[] = "make_clov";

  f_mn = (su3_matrix *)malloc(sizeof(su3_matrix)*sites_on_node);
  if(f_mn == NULL){
    printf("%s(%d): Can't malloc f_mn\n",myname,this_node);
    terminate(1);
  }



/* The clover term
*                                         ( X | 0 )
* sum_{mu<nu} sigma_{mu,nu} F_{mu,nu} = ( ----- )
*                                         ( 0 | Y )
*
* is block diagonal in the Weyl basis used, and hermitian. Together with the 1|
* we will store it in a lower complex triangular (block) matrix (without
* diagonal) and a real diagonal. The blocks go over Dirac indices 0,1 and 2,3.
*
*   With
*
*   sigma_{mu,nu} = i/2 [gamma_{mu}, gamma_{nu}]
*
*   and with our gamma matrices (see libraries/mb_gamma.c)
*
*   1234 = XYZT    and   Pauli sigma_1, sigma_2, sigma_3
*
*   sigma_23 =  diag( sigma_1,  sigma_1)
*   sigma_31 =  diag(-sigma_2, -sigma_2)
*   sigma_12 =  diag( sigma_3,  sigma_3)
*   sigma_14 =  diag(-sigma_1,  sigma_1)
*   sigma_24 =  diag( sigma_2, -sigma_2)
*   sigma_34 =  diag(-sigma_3,  sigma_3)
*
*
* Here X = sigma_1 (F_23 - F_14) + sigma_2 (F_13 + F_24) + sigma_3 (F_12 - F_34)
* and  Y = sigma_1 (F_23 + F_14) + sigma_2 (F_13 - F_24) + sigma_3 (F_12 + F_34)
*
* This has to be multiplied by the "clover coefficient" and subtracted from 1|.
*
* Note that F_mn = f_mn/i = -i f_mn!
*
* (The indices above start from 1, in the program they start from 0,
* i.e. subtract -1 from the above) */

/* The clover term R is a hermitian 12 x 12 matrix in color and Dirac spin.
   If we use ab for color and ij for spin, with a = (0,1,2) and
   i = (0,1,2,3), then we store the matrix R(ai,bj) as follows with
   the color index varying most rapidly:

   Upper left block, showing blocks of SU(3) matrices

   [ 1 + ci(f_01 - f_23)                ci(f_12 - f_03) + c(f_02 + f_13) ]
   [ ci(f_12 - f_03) - c(f_02 + f_13)   1 - ci(f_01 - f_23)                ]

   Lower right block, showing blocks of SU(3) matrices

   [ 1 + ci(f_01 + f_23)                 ci(f_12 + f_03) + c(f_02 - f_13)  ]
   [ ci(f_12 + f_03) - c(f_02 - f_13)    1 - ci(f_01 + f_23)               ]

   This hermitian matrix is decomposed into a lower triangular matrix
   and its diagonals.  The diagonals are real.

   clov_diag[site].di[0][0-5] holds the diagonals from the upper left
   clov_diag[site].di[1][0-5] holds the diagonals from the lower right

   clov_diag[site].di[0][a]   = R(a0,a0) = 1 - c(Im f_01 - Im f_23)(a,a)
   clov_diag[site].di[0][a+3] = R(a1,a1) = 1 + c(Im f_01 - Im f_23)(a,a)
   clov_diag[site].di[1][a]   = R(a2,a2) = 1 - c(Im f_01 + Im f_23)(a,a)
   clov_diag[site].di[1][a+3] = R(a3,a3) = 1 + c(Im f_01 + Im f_23)(a,a)

   clov[site].tr[0][jk]:
      Holds the lower triangular elements of the upper left block in
      the pattern

        jk = ...

   [          |        ]      [                      |                     ]
   [  0       |        ]      [  1 + ci(f_01 - f_23) | ci(f_12 - f_03)     ]
   [  1  2    |        ]      [                      |   + c(f_02 + f_13)  ]
   ---------------------   =  ----------------------------------------------
   [  3  4  5 |        ]      [			     |		           ]
   [  6  7  8 |  9     ]      [ ci(f_12 - f_03)      | 1 - ci(f_01 - f_23) ]
   [ 10 11 12 | 13 14  ]      [    - c(f_02 + f_13)  |		           ]

   clov[site].tr[1][jk]:
      Holds the lower triangular elements of the lower right block in
      the pattern

        jk = ...

   [          |        ]      [                      |                     ]
   [  0       |        ]      [  1 + ci(f_01 + f_23) | ci(f_12 + f_03)     ]
   [  1  2    |        ]      [                      |   + c(f_02 - f_13)  ]
   ---------------------   =  ----------------------------------------------
   [  3  4  5 |        ]      [			     |		           ]
   [  6  7  8 |  9     ]      [ ci(f_12 + f_03)      | 1 - ci(f_01 + f_23) ]
   [ 10 11 12 | 13 14  ]      [    - c(f_02 - f_13)  |		           ]

   */

    f_mu_nu(f_mn, 0, 1);
    FORALLSITESDOMAIN(i,s){
        scalar_mult_su3_matrix( f_mn+i, Clov_c,
	    (f_mn+i) );
    }
    jk=0;
    for(j=0;j<DIMF;j++){
	FORALLSITESDOMAIN(i,s){
	    clov_diag[i].di[0][j] = 1.0 -
		(f_mn+i)->e[j][j].imag;
	    clov_diag[i].di[0][j+DIMF] = 1.0 +
		(f_mn+i)->e[j][j].imag;
	    clov_diag[i].di[1][j] = 1.0 -
		(f_mn+i)->e[j][j].imag;
	    clov_diag[i].di[1][j+DIMF] = 1.0 +
		(f_mn+i)->e[j][j].imag;
	}
	jk2=(j+DIMF)*(j+DIMF-1)/2+DIMF;
	for(k=0;k<j;k++){
	    FORALLSITESDOMAIN(i,s){
		TIMESPLUSI( (f_mn+i)->e[j][k],
		    clov[i].tr[0][jk]);
		TIMESMINUSI( (f_mn+i)->e[j][k],
		    clov[i].tr[0][jk2]);
		TIMESPLUSI( (f_mn+i)->e[j][k],
		    clov[i].tr[1][jk]);
		TIMESMINUSI( (f_mn+i)->e[j][k],
		    clov[i].tr[1][jk2]);
	    }
	    jk++;
	    jk2++;
	}
    }

    f_mu_nu(f_mn, 2, 3);
    FORALLSITESDOMAIN(i,s){
        scalar_mult_su3_matrix( (f_mn+i), Clov_c,
	    (f_mn+i) );
    }
    jk=0;
    for(j=0;j<DIMF;j++){
	FORALLSITESDOMAIN(i,s){
	    clov_diag[i].di[0][j] +=
		(f_mn+i)->e[j][j].imag;
	    clov_diag[i].di[0][j+DIMF] -=
		(f_mn+i)->e[j][j].imag;
	    clov_diag[i].di[1][j] -=
		(f_mn+i)->e[j][j].imag;
	    clov_diag[i].di[1][j+DIMF] +=
		(f_mn+i)->e[j][j].imag;
	}
	jk2=(j+DIMF)*(j+DIMF-1)/2+DIMF;
	for(k=0;k<j;k++){
	    FORALLSITESDOMAIN(i,s){
		TIMESMINUSI( (f_mn+i)->e[j][k], ctmp);
		CADD( clov[i].tr[0][jk], ctmp, clov[i].tr[0][jk]);
		CSUB( clov[i].tr[0][jk2], ctmp, clov[i].tr[0][jk2]);
		CSUB( clov[i].tr[1][jk], ctmp, clov[i].tr[1][jk]);
		CADD( clov[i].tr[1][jk2], ctmp, clov[i].tr[1][jk2]);
	    }
	    jk++;
	    jk2++;
	}
    }

    f_mu_nu(f_mn, 1, 2);
    FORALLSITESDOMAIN(i,s){
        scalar_mult_su3_matrix( (f_mn+i), Clov_c,
	    (f_mn+i) );
    }
    for(j=0;j<DIMF;j++){
	jk=(j+DIMF)*(j+DIMF-1)/2;
	for(k=0;k<DIMF;k++){
	    FORALLSITESDOMAIN(i,s){
		TIMESPLUSI( (f_mn+i)->e[j][k],
		    clov[i].tr[0][jk]);
		TIMESPLUSI( (f_mn+i)->e[j][k],
		    clov[i].tr[1][jk]);
	    }
	    jk++;
	}
    }

    f_mu_nu(f_mn, 0, 3);
    FORALLSITESDOMAIN(i,s){
        scalar_mult_su3_matrix( (f_mn+i), Clov_c,
	    (f_mn+i) );
    }
    for(j=0;j<DIMF;j++){
	jk=(j+DIMF)*(j+DIMF-1)/2;
	for(k=0;k<DIMF;k++){
	    FORALLSITESDOMAIN(i,s){
		TIMESMINUSI( (f_mn+i)->e[j][k], ctmp);
		CADD( clov[i].tr[0][jk], ctmp, clov[i].tr[0][jk]);
		CSUB( clov[i].tr[1][jk], ctmp, clov[i].tr[1][jk]);
	    }
	    jk++;
	}
    }

    f_mu_nu(f_mn, 0, 2);
    FORALLSITESDOMAIN(i,s){
        scalar_mult_su3_matrix( (f_mn+i), Clov_c,
	    (f_mn+i) );
    }
    for(j=0;j<DIMF;j++){
	jk=(j+DIMF)*(j+DIMF-1)/2;
	for(k=0;k<DIMF;k++){
	    FORALLSITESDOMAIN(i,s){
		CSUB( clov[i].tr[0][jk],
		    (f_mn+i)->e[j][k], clov[i].tr[0][jk]);
		CSUB( clov[i].tr[1][jk],
		    (f_mn+i)->e[j][k], clov[i].tr[1][jk]);
	    }
	    jk++;
	}
    }

    f_mu_nu(f_mn, 1, 3);
    FORALLSITESDOMAIN(i,s){
        scalar_mult_su3_matrix( (f_mn+i), Clov_c,
	    (f_mn+i) );
    }
    for(j=0;j<DIMF;j++){
	jk=(j+DIMF)*(j+DIMF-1)/2;
	for(k=0;k<DIMF;k++){
	    FORALLSITESDOMAIN(i,s){
		CSUB( clov[i].tr[0][jk],
		    (f_mn+i)->e[j][k], clov[i].tr[0][jk]);
		CADD( clov[i].tr[1][jk],
		    (f_mn+i)->e[j][k], clov[i].tr[1][jk]);
	    }
	    jk++;
	}
    }

    free(f_mn);

} /* compute_clov */

/******************* compute_clovinv ***********************************/

/* MIMD version 7 */

/* Modifications
   added parity selection CD
   version of 1/3/95 by UMH
*/

/* Invert the "clover" term on a specified sublattice and return tr(log(A)) */
/* (tr log needed for clover_dynamical) */

double compute_clovinv(clover *my_clov, int parity) {

  triangular *clov = my_clov->clov;
  diagonal *clov_diag = my_clov->clov_diag;
  register int i;
  register site *s;

  register int b,j,k,kj,l,lk,lj,jl;
  Real f_diag[2*DIMF];
  complex v1[2*DIMF],ctmp,sum;
  double trlogA;

    /* Take the inverse on the odd sublattice for each of the 2 blocks */
    trlogA = (double)0.0;
    FORSOMEPARITYDOMAIN(i,s,parity)for(b=0;b<2;b++){

	/* Cholesky decompose. This can be done in place, except that
	   we will need a temporary diagonal.
	   Uses algorithm 4.2.2 of "Matrix Computations" by Gene H. Golub
	   and Charles F. Van Loan (John Hopkins, 1989) */
	for(j=0;j<2*DIMF;j++){

	    clov_diag[i].di[b][j] = sqrt((double)(clov_diag[i].di[b][j]));
	    f_diag[j] = 1.0 / (clov_diag[i].di[b][j]);
	    for(k=j+1;k<2*DIMF;k++){
		kj=k*(k-1)/2+j;
		CMULREAL(clov[i].tr[b][kj], f_diag[j], clov[i].tr[b][kj]);
	    }

	    for(k=j+1;k<2*DIMF;k++){
		kj=k*(k-1)/2+j;
		CMUL_J(clov[i].tr[b][kj], clov[i].tr[b][kj], ctmp);
		clov_diag[i].di[b][k] -= ctmp.real;
		for(l=k+1;l<2*DIMF;l++){
		    lj=l*(l-1)/2+j;
		    lk=l*(l-1)/2+k;
		    CMUL_J(clov[i].tr[b][lj], clov[i].tr[b][kj], ctmp);
		    CSUB(clov[i].tr[b][lk], ctmp, clov[i].tr[b][lk]);
		}
	    }
	}

	/* Accumulate trlogA */
	for(j=0;j<2*DIMF;j++)
	    trlogA += (double)2.0 * log((double)(clov_diag[i].di[b][j]));

	/* Now use forward and backward substitution to construct inverse */
	for(k=0;k<2*DIMF;k++){
	    for(l=0;l<k;l++) v1[l] = cmplx(0.0, 0.0);

	    /* Forward substitute */
	    v1[k] = cmplx(f_diag[k], 0.0);
	    for(l=k+1;l<2*DIMF;l++){
		sum = cmplx(0.0, 0.0);
		for(j=k;j<l;j++){
		    lj=l*(l-1)/2+j;
		    CMUL(clov[i].tr[b][lj], v1[j], ctmp);
		    CSUB(sum, ctmp, sum);
		}
		CMULREAL(sum, f_diag[l], v1[l]);
	    }

	    /* Backward substitute */
	    CMULREAL(v1[2*DIMF-1], f_diag[2*DIMF-1], v1[2*DIMF-1]);
	    for(l=2*DIMF-2;l>=k;l--){
		sum = v1[l];
		for(j=l+1;j<2*DIMF;j++){
		    jl=j*(j-1)/2+l;
		    CMULJ_(clov[i].tr[b][jl], v1[j], ctmp);
		    CSUB(sum, ctmp, sum);
		}
		CMULREAL(sum, f_diag[l], v1[l]);
	    }

	    /* Overwrite column k */
	    clov_diag[i].di[b][k] = v1[k].real;
	    for(l=k+1;l<2*DIMF;l++){
		lk=l*(l-1)/2+k;
		clov[i].tr[b][lk] = v1[l];
	    }
	}
    } END_LOOP

    g_doublesum( &trlogA );
    return(trlogA);

} /* compute_clovinv */

/******************* mult_this_ldu_site *********************************/

/* version of 12/29/94 by UMH */
/* CD 06/02/98 loop rearrangements, prefetch and macros as suggested by DT */

/* Multiply a Wilson vector (spin,color) with a block-diagonal hermition
   matrix stored as a complex lower triangular matrix (without diagonal)
   and a real diagonal. The blocks are spin indices 0,1 and 2,3. */

/* To simplify the code with the above block structure, introduce
   the somewhat dirty structure, equivalent to a wilson_vector: */
typedef struct { complex b[2][2*DIMF]; } wilson_block_vector;

/* c += a*b */
#define CMULSUM(a,b,c) { (c).real += (a).real*(b).real - (a).imag*(b).imag; \
				(c).imag += (a).real*(b).imag + (a).imag*(b).real; }
/*  c += conj(a)*b */
#define CMULJ_SUM(a,b,c) { (c).real += (a).real*(b).real + (a).imag*(b).imag; \
				(c).imag += (a).real*(b).imag - (a).imag*(b).real; }

void mult_this_ldu_site(
  clover *my_clov,
  field_offset src,   /* type wilson_vector RECAST AS wilson_block_vector */
  field_offset dest,  /* type wilson_vector RECAST AS wilson_block_vector */
  int parity
  )
{
  triangular *clov = my_clov->clov;
  diagonal *clov_diag = my_clov->clov_diag;

  register int i;
  register site *s;
  register int b,j,k,jk,kj;
  register complex ctmp;

  FORSOMEPARITYDOMAIN(i,s,parity){
    if( i < loopend-FETCH_UP ){
      prefetch_W( (wilson_vector *)(F_PT((s+FETCH_UP),src)) );
      prefetch_W( (wilson_vector *)(F_PT((s+FETCH_UP),dest)) );
    }
    for(b=0;b<2;b++){
      for(j=0;j<2*DIMF;j++){

	/* diagonal part */
	CMULREAL(((wilson_block_vector *)F_PT(s,src))->b[b][j],
		 clov_diag[i].di[b][j],
		 ctmp );

	/* lower triangular part */
	jk=j*(j-1)/2;
	for(k=0;k<j;k++){
	  CMULSUM(clov[i].tr[b][jk],
		  ((wilson_block_vector *)F_PT(s,src))->b[b][k],
		  ctmp );
	  jk++;
	}

	/* upper triangular part */
	for(k=j+1;k<2*DIMF;k++){
	  kj=k*(k-1)/2+j;
	  CMULJ_SUM(clov[i].tr[b][kj],
		    ((wilson_block_vector *)F_PT(s,src))->b[b][k],
		    ctmp );
	}
	((wilson_block_vector *)F_PT(s,dest))->b[b][j] = ctmp;
      }  /* j loop */
    } /* b loop */
  } END_LOOP /* site loop */

} /* mult_this_ldu_site */


/******************* mult_this_ldu_field **********************************/

void mult_this_ldu_field(
  clover *my_clov,
  wilson_vector *src,  /* RECAST as wilson_block_vector */
  wilson_vector *dest, /* RECAST as wilson_block_vector */
  int parity
  )
{
  triangular *clov = my_clov->clov;
  diagonal *clov_diag = my_clov->clov_diag;
  register int i;
  register site *s;
  register int b,j,k,jk,kj;
  register complex ctmp;
  wilson_block_vector *srcb = (wilson_block_vector *)src;
  wilson_block_vector *destb = (wilson_block_vector *)dest;

  FORSOMEPARITYDOMAIN(i,s,parity){
    if( i < loopend-FETCH_UP ){
      prefetch_W( &srcb[i+FETCH_UP] );
      prefetch_W( &destb[i+FETCH_UP] );
    }
    for(b=0;b<2;b++){
      for(j=0;j<2*DIMF;j++){

	/* diagonal part */
	CMULREAL( srcb[i].b[b][j], clov_diag[i].di[b][j], ctmp );

	/* lower triangular part */
	jk=j*(j-1)/2;
	for(k=0;k<j;k++){
	  CMULSUM(clov[i].tr[b][jk],
		  srcb[i].b[b][k],
		  ctmp );
	  jk++;
	}

	/* upper triangular part */
	for(k=j+1;k<2*DIMF;k++){
	  kj=k*(k-1)/2+j;
	  CMULJ_SUM(clov[i].tr[b][kj], srcb[i].b[b][k], ctmp );
	}
	destb[i].b[b][j] = ctmp;
      }  /* j loop */
    } /* b loop */
  } END_LOOP /* site loop */

} /* mult_this_ldu_field */

/******** tr_sigma_ldu_site *************/
/* MIMD version 7, May 95, MBW */
/* Used in clover_diagonal */
/* Left multiplies a color-dirac matrix stored in "ldu" form
by sigma_mu_nu = gamma_mu gamma_nu where

sigma(XUP,YUP)          sigma(XUP,ZUP)          sigma(XUP,TUP)
        -i  0  0  0              0 -1  0  0              0  i  0  0
         0  i  0  0              1  0  0  0              i  0  0  0
         0  0 -i  0              0  0  0 -1              0  0  0 -i
         0  0  0  i              0  0  1  0              0  0 -i  0

sigma(YUP,ZUP)          sigma(YUP,TUP)          sigma(ZUP,TUP)
         0 -i  0  0              0 -1  0  0              i  0  0  0
        -i  0  0  0              1  0  0  0              0 -i  0  0
         0  0  0 -i              0  0  0  1              0  0 -i  0
         0  0 -i  0              0  0 -1  0              0  0  0  i

and sigma(nu,mu) = -sigma(mu,nu), and sums over the dirac indices.

*/


#define XUP 0
#define YUP 1
#define ZUP 2
#define TUP 3

void tr_sigma_ldu_mu_nu_site( field_offset mat, int mu, int nu )
/* triang, diag are input & contain the color-dirac matrix
   mat is output: the resulting su3_matrix  */
{
  triangular *clov = global_clov->clov;
  diagonal *clov_diag = global_clov->clov_diag;
  register site *s;
  register int mm = 0, nn = 0;  /* dummy directions */
  register Real pm = 0;
  register int i,j,k,jk,jk2,kj;
  Real rtmp;
  complex ctmp,ctmp1;
  char myname[] = "tr_sigma_ldu_mu_nu_site";

    /* take care of the case mu > nu by flipping them and mult. by -1 */
    if( mu < nu ) {
	mm = mu; nn = nu; pm = 1.0;
    }
    else if( mu > nu) {
	mm = nu; nn = mu; pm = -1.0;
    }
    else
	printf("BAD CALL by %s: mu=%d, nu=%d\n",myname,mu,nu);

switch(mm) {
    case(XUP):
    switch(nn) {
	case(YUP):
	FORODDSITESDOMAIN(i,s) {
	    /* diagonal part */
	    for(j=0;j<DIMF;j++) {
		rtmp = clov_diag[i].di[0][j+DIMF]
		    - clov_diag[i].di[0][j];
		rtmp -= clov_diag[i].di[1][j];
		rtmp += clov_diag[i].di[1][j+DIMF];
		((su3_matrix *)F_PT(s,mat))->e[j][j].real = 0.0;
		((su3_matrix *)F_PT(s,mat))->e[j][j].imag = rtmp;
	    }
	    /* triangular part */
	    jk = 0;
	    for(j=1;j<DIMF;j++) {
		jk2 = (j+DIMF)*(j+DIMF-1)/2 + DIMF;
		for(k=0;k<j;k++) {
		    CSUB( clov[i].tr[0][jk2],
			clov[i].tr[0][jk], ctmp );
		    CSUB( ctmp, clov[i].tr[1][jk],
			ctmp );
		    CADD( ctmp, clov[i].tr[1][jk2],
			ctmp );
		    TIMESPLUSI( ctmp, ((su3_matrix *)F_PT(s,mat))->e[j][k] );
		    CONJG( ctmp, ctmp);
		    TIMESPLUSI( ctmp, ((su3_matrix *)F_PT(s,mat))->e[k][j] );
		    jk++; jk2++;
		}
	    }
	}
	break;
	case(ZUP):
	FORODDSITESDOMAIN(i,s) {
	    for(j=0;j<DIMF;j++) {
		jk = (j+DIMF)*(j+DIMF-1)/2;
		for(k=0;k<DIMF;k++) {
		    kj = (k+DIMF)*(k+DIMF-1)/2 + j;
		    CONJG( clov[i].tr[0][kj], ctmp);
		    CSUB( ctmp, clov[i].tr[0][jk],				ctmp );
		    CONJG( clov[i].tr[1][kj], ctmp1);
		    CADD( ctmp, ctmp1, ctmp);
		    CSUB( ctmp, clov[i].tr[1][jk],
			((su3_matrix *)F_PT(s,mat))->e[j][k] );
		    jk++;
		}
	    }
	}
	break;
	case(TUP):
	FORODDSITESDOMAIN(i,s) {
	    for(j=0;j<DIMF;j++) {
		jk = (j+DIMF)*(j+DIMF-1)/2;
		for(k=0;k<DIMF;k++) {
		    kj = (k+DIMF)*(k+DIMF-1)/2 + j;
		    CONJG( clov[i].tr[0][kj], ctmp);
		    CADD( ctmp, clov[i].tr[0][jk],
			ctmp );
		    CONJG( clov[i].tr[1][kj], ctmp1);
		    CSUB( ctmp, ctmp1, ctmp);
		    CSUB( ctmp, clov[i].tr[1][jk],
			ctmp );
		    TIMESPLUSI( ctmp, ((su3_matrix *)F_PT(s,mat))->e[j][k] );
		    jk++;
		}
	    }
	}
	break;
	default:
	   printf("BAD CALL in %s: mu=%d, nu=%d\n",myname,mu,nu);
    }
    break;
    case(YUP):
    switch(nn) {
	case(ZUP):
	FORODDSITESDOMAIN(i,s) {
	    for(j=0;j<DIMF;j++) {
		jk = (j+DIMF)*(j+DIMF-1)/2;
		for(k=0;k<DIMF;k++) {
		    kj = (k+DIMF)*(k+DIMF-1)/2 + j;
		    CONJG( clov[i].tr[0][kj], ctmp);
		    CADD( ctmp, clov[i].tr[0][jk],
			ctmp );
		    CONJG( clov[i].tr[1][kj], ctmp1);
		    CADD( ctmp, ctmp1, ctmp);
		    CADD( ctmp, clov[i].tr[1][jk],
			ctmp );
		    TIMESMINUSI( ctmp, ((su3_matrix *)F_PT(s,mat))->e[j][k] );
		    jk++;
		}
	    }
	}
	break;
	case(TUP):
	FORODDSITESDOMAIN(i,s) {
	    for(j=0;j<DIMF;j++) {
		jk = (j+DIMF)*(j+DIMF-1)/2;
		for(k=0;k<DIMF;k++) {
		    kj = (k+DIMF)*(k+DIMF-1)/2 + j;
		    CONJG( clov[i].tr[0][kj], ctmp);
		    CSUB( ctmp, clov[i].tr[0][jk],
			ctmp );
		    CONJG( clov[i].tr[1][kj], ctmp1);
		    CSUB( ctmp, ctmp1, ctmp);
		    CADD( ctmp, clov[i].tr[1][jk],
			((su3_matrix *)F_PT(s,mat))->e[j][k] );
		    jk++;
		}
	    }
	}
	break;
	default:
	   printf("BAD CALL in %s: mu=%d, nu=%d\n",myname,mu,nu);
    }
    break;
    case(ZUP):
    switch(nn) {
	case(TUP):
	FORODDSITESDOMAIN(i,s) {
	    /* diagonal part */
	    for(j=0;j<DIMF;j++) {
		rtmp = clov_diag[i].di[0][j]
		    - clov_diag[i].di[0][j+DIMF];
		rtmp -= clov_diag[i].di[1][j];
		rtmp += clov_diag[i].di[1][j+DIMF];
		((su3_matrix *)F_PT(s,mat))->e[j][j].real = 0.0;
		((su3_matrix *)F_PT(s,mat))->e[j][j].imag = rtmp;
	    }
	    /* triangular part */
	    jk = 0;
	    for(j=1;j<DIMF;j++) {
		jk2 = (j+DIMF)*(j+DIMF-1)/2 + DIMF;
		for(k=0;k<j;k++) {
		    CSUB( clov[i].tr[0][jk],
			clov[i].tr[0][jk2], ctmp );
		    CSUB( ctmp, clov[i].tr[1][jk],
			ctmp );
		    CADD( ctmp, clov[i].tr[1][jk2],
			ctmp );
		    TIMESPLUSI( ctmp, ((su3_matrix *)F_PT(s,mat))->e[j][k] );
		    CONJG( ctmp, ctmp);
		    TIMESPLUSI( ctmp, ((su3_matrix *)F_PT(s,mat))->e[k][j] );
		    jk++; jk2++;
		}
	    }
	}
	break;
	default:
	   printf("BAD CALL in %s: mu=%d, nu=%d\n",myname,mu,nu);
    }
    break;
    default:
	printf("BAD CALL in %s: mu=%d, nu=%d\n",myname,mu,nu);
}

/* multiply by pm = +/- 1  */
FORODDSITESDOMAIN(i,s) {
    for(j=0;j<DIMF;j++) for(k=0;k<DIMF;k++)
	CMULREAL( ((su3_matrix *)F_PT(s,mat))->e[j][k], pm,
	    ((su3_matrix *)F_PT(s,mat))->e[j][k] );
}

} /* tr_sigma_ldu */


/******************* free_this_clov *********************************/

void free_this_clov(clover *my_clov)
{

  free(my_clov->clov);
  free(my_clov->clov_diag);
  free(my_clov);

} /* free_this_clov */


/******************* make_clov ***************************************/
/* For making only the global clover term */

void make_clov(Real Clov_c)
{
  global_clov = create_clov();
  if(global_clov == NULL)
    terminate(1);

  compute_clov(global_clov, Clov_c);
}

/******************* make_clovinv ***********************************/
/* For inverting only the global clover term */

double make_clovinv(int parity)
{
  return compute_clovinv(global_clov, parity);
}


/******************* mult_ldu_site ***********************************/
/* For multiplying only by the global clover term */
void mult_ldu_site(
  field_offset src,   /* type wilson_vector RECAST AS wilson_block_vector */
  field_offset dest,  /* type wilson_vector RECAST AS wilson_block_vector */
  int parity
  )
{
  mult_this_ldu_site(global_clov, src, dest, parity);
}

/******************* mult_ldu_field ***********************************/
/* For multiplying only by the global clover term */
void mult_ldu_field(
  wilson_vector *src,
  wilson_vector *dest,
  int parity
  )
{
  mult_this_ldu_field(global_clov, src, dest, parity);
}

/******************* free_clov ***************************************/

void free_clov()
{
  free_this_clov(global_clov);

} /* free_clov */

