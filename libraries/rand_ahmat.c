/******************  rand_ahmat.c  (in su3.a) ***************************
*									*
* void random_anti_hermitian( anti_hermitmat *mat_antihermit, passthru *prn_pt)*
* Creates gaussian random anti-hermitian matrices			*
* Normalization is < |m01|^2 > = 1, or < m01.real*m01.real > = 1/2	*
* The argument "prn_pt" is a pointer to be passed to gaussian_rand_no() *
* RS6000 may choke on void *						*
*/
#include "../include/config.h"
#include <math.h>
#include "../include/complex.h"
#include "../include/su3.h"

void random_anti_hermitian( anti_hermitmat *mat_antihermit, double_prn *prn_pt) {
Real r3;
#if (NCOL>2)
Real r8,sqrt_third;
#if (NCOL>3)
Real r15,sqrt_sixth;
#endif
#endif

	/* off-diagonal elements */
        r3=gaussian_rand_no(prn_pt);
	mat_antihermit->m01.real=gaussian_rand_no(prn_pt);
	mat_antihermit->m01.imag=gaussian_rand_no(prn_pt);
#if (NCOL>2)
	sqrt_third = sqrt( (double)(1.0/3.0) );
	r8=gaussian_rand_no(prn_pt);
	mat_antihermit->m02.real=gaussian_rand_no(prn_pt);
	mat_antihermit->m12.real=gaussian_rand_no(prn_pt);
	mat_antihermit->m02.imag=gaussian_rand_no(prn_pt);
	mat_antihermit->m12.imag=gaussian_rand_no(prn_pt);
#if (NCOL>3)
	sqrt_sixth = sqrt( (double)(1.0/6.0) );
	r15=gaussian_rand_no(prn_pt);
	mat_antihermit->m03.real=gaussian_rand_no(prn_pt);
	mat_antihermit->m13.real=gaussian_rand_no(prn_pt);
	mat_antihermit->m23.real=gaussian_rand_no(prn_pt);
	mat_antihermit->m03.imag=gaussian_rand_no(prn_pt);
	mat_antihermit->m13.imag=gaussian_rand_no(prn_pt);
	mat_antihermit->m23.imag=gaussian_rand_no(prn_pt);
#endif
#endif
	/* diagonal elements */
#if (NCOL==2)
	mat_antihermit->m00im=  r3;
	mat_antihermit->m11im= -r3;
#endif
#if (NCOL==3)
	mat_antihermit->m00im=   r3+sqrt_third*r8;
	mat_antihermit->m11im=  -r3+sqrt_third*r8;
	mat_antihermit->m22im= -2.0*sqrt_third*r8;
#endif
#if (NCOL==4)
	mat_antihermit->m00im=   r3+sqrt_third*r8+sqrt_sixth*r15;
	mat_antihermit->m11im=  -r3+sqrt_third*r8+sqrt_sixth*r15;
	mat_antihermit->m22im= -2.0*sqrt_third*r8+sqrt_sixth*r15;
	mat_antihermit->m33im=               -3.0*sqrt_sixth*r15;
#endif

}/*random_anti_hermitian_*/
