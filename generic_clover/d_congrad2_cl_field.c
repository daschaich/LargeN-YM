/******* conjugate gradient for clover fermions ****/
/* Clover fermions */
/* For clover_dynamical/update.c.  Solves M_adjoint*M psi = chi */

/* if "LU" is defined use the LU preconditioned fermion matrix, where
  the fermion spinors live on even sites only.  In other words, if
  Dslash_oe is the dslash operator with its source on even sites and
  its result on odd sites, etc.:

  without LU:
   M = A - kappa*( Dslash_eo + DSLASH_oe )
  with LU:
   M = A_e - kappa^2 * Dslash_eo * (A_o)^{-1} * Dslash_oe
*/

/* adapted to "field" instead of "site" */

#ifdef LU
#define FORMYSITESDOMAIN FOREVENSITESDOMAIN
#else
#define FORMYSITESDOMAIN FORALLSITESDOMAIN
#endif

/* This version looks at the initial vector every "niter" passes */
/* The source vector is in "chi", and the initial guess and answer
   in "psi".  "r" is the residual vector, and "p" and "mp" are
   working vectors for the conjugate gradient.
   niter = maximum number of iterations.
   rsqmin = desired rsq, quit when we reach rsq = rsqmin*source_norm.
*/

#include "generic_clover_includes.h"

int congrad_cl(int niter,Real rsqmin,Real *final_rsq_ptr) {
register int i;
register site *s;
int iteration;	/* counter for iterations */
double source_norm;
double rsqstop;
Real a,b;
double rsq,oldrsq,pkp;	/* Sugar's a,b,resid**2,previous resid*2 */
				/* pkp = cg_p.K.cg_p */
void dslash_w_field_special();
msg_tag *tag[8],*tag2[8];
#ifdef LU
Real KAP = -kappa*kappa;
#else
Real KAP = -kappa;
#endif

/* create space for fields, copy source and initial guess */

wilson_vector *psi, *chi, *tmp, *mp, *p, *r;
    FIELD_ALLOC(psi,wilson_vector)
    FIELD_ALLOC(chi,wilson_vector)
    FIELD_ALLOC(tmp,wilson_vector)
    FIELD_ALLOC(mp,wilson_vector)
    FIELD_ALLOC(p,wilson_vector)
    FIELD_ALLOC(r,wilson_vector)

    FORALLSITES(i,s){
        psi[i] = s->psi;
        chi[i] = s->chi;
    }

double dtime;
dtime= -dclock();

	iteration=0;
start:
	/* mp <-  M_adjoint*M*psi
	   r,p <- chi - mp
	   rsq = |r|^2
	   source_norm = |chi|^2
	*/
	rsq = source_norm = 0.0;
#ifdef LU
	mult_ldu_field(psi, tmp, EVEN);
        dslash_w_field_special(psi, mp, PLUS, ODD, tag, 0);
	mult_ldu_field(mp, tmp, ODD);
        dslash_w_field_special(tmp, mp, PLUS, EVEN, tag2, 0);
        FOREVENSITESDOMAIN(i,s){
            scalar_mult_add_wvec( tmp+i, mp+i, KAP, mp+i);
        }
	mult_ldu_field(mp, tmp, EVEN);
        dslash_w_field_special(mp, mp, MINUS, ODD, tag, 1);
	mult_ldu_field(mp, tmp, ODD);
        dslash_w_field_special(tmp, mp, MINUS, EVEN, tag2, 1);
        FOREVENSITESDOMAIN(i,s){
            scalar_mult_add_wvec( tmp+i, mp+i, KAP, mp+i );
	    sub_wilson_vector( chi+i, mp+i, r+i );
	    p[i] = r[i];
	    source_norm += (double)magsq_wvec( chi+i );
	    rsq += (double)magsq_wvec( r+i );
        }
#else
	mult_ldu_field(psi, tmp, EVENANDODD);
	dslash_w_field_special(psi, mp, PLUS, EVENANDODD, tag, 0);
	FORALLSITESDOMAIN(i,s){
	    scalar_mult_add_wvec( tmp+i, mp+i, KAP, mp+i );
	}
	mult_ldu_field(mp, tmp, EVENANDODD);
	dslash_w_field_special(mp, mp, MINUS, EVENANDODD, tag, 1);
	FORALLSITESDOMAIN(i,s){
	    scalar_mult_add_wvec( tmp+i, mp+i, KAP, mp+i );
	    sub_wilson_vector( chi+i, mp+i, r+i );
	    p[i] = r[i];
	    source_norm += (double)magsq_wvec( chi+i );
	    rsq += (double)magsq_wvec( r+i );
	}
#endif
	g_doublesum( &source_norm );
	g_doublesum( &rsq );
        iteration++ ;	/* iteration counts number of multiplications
			   by M_adjoint*M */
	total_iters++;
/**if(this_node==0)printf("congrad2: source_norm = %e\n",source_norm);
if(this_node==0)printf("congrad2: iter %d, rsq %e, pkp %e, a %e\n",
iteration,(double)rsq,(double)pkp,(double)a );**/

	rsqstop = rsqmin * source_norm;
	if( rsq <= rsqstop ){
	    *final_rsq_ptr= (Real)rsq;
	    for( i=XUP; i <= TUP; i++) {
		cleanup_gather(tag[i]);
		cleanup_gather(tag[OPP_DIR(i)]);
#ifdef LU
		cleanup_gather(tag2[i]);
		cleanup_gather(tag2[OPP_DIR(i)]);
#endif
	    }

            cleanup_dslash_temps();
            cleanup_tmp_links();

            FORALLSITES(i,s){
                s->psi = psi[i];
            }
            free(psi); free(chi); free(tmp);
            free(mp); free(p); free(r);

	    return (iteration);
	}

    /* main loop - do until convergence or time to restart */
	/*
	   oldrsq <- rsq
	   mp <- M_adjoint*M*p
	   pkp <- p.M_adjoint*M.p
	   a <- rsq/pkp
	   psi <- psi + a*p
	   r <- r - a*mp
	   rsq <- |r|^2
	   b <- rsq/oldrsq
	   p <- r + b*p
	*/
    do{
	oldrsq = rsq;
	pkp = 0.0;
#ifdef LU
	mult_ldu_field(p, tmp, EVEN);
        dslash_w_field_special(p, mp, PLUS, ODD, tag, 1);
	mult_ldu_field(mp, tmp, ODD);
        dslash_w_field_special(tmp, mp, PLUS, EVEN, tag2, 1);
        FOREVENSITESDOMAIN(i,s){
            scalar_mult_add_wvec( tmp+i, mp+i, KAP, mp+i );
        }
	mult_ldu_field(mp, tmp, EVEN);
        dslash_w_field_special(mp, mp, MINUS, ODD, tag, 1);
	mult_ldu_field(mp, tmp, ODD);
        dslash_w_field_special(tmp, mp, MINUS, EVEN, tag2, 1);
        FOREVENSITESDOMAIN(i,s){
            scalar_mult_add_wvec( tmp+i, mp+i, KAP, mp+i );
            pkp += (double)wvec_rdot( p+i, mp+i );
        }
#else
	mult_ldu_field(p, tmp, EVENANDODD);
	dslash_w_field_special(p, mp, PLUS, EVENANDODD, tag, 1);
	FORALLSITESDOMAIN(i,s){
	    scalar_mult_add_wvec( tmp+i, mp+i, KAP, mp+i );
	}
	mult_ldu_field(mp, tmp, EVENANDODD);
	dslash_w_field_special(mp, mp, MINUS, EVENANDODD, tag, 1);
	FORALLSITESDOMAIN(i,s){
	    scalar_mult_add_wvec( tmp+i, mp+i, KAP, mp+i );
            pkp += (double)wvec_rdot( p+i, mp+i );
	}
#endif
	g_doublesum( &pkp );
	iteration++;
	total_iters++;

	a = (Real)(rsq/pkp);
	rsq = 0.0;
	FORMYSITESDOMAIN(i,s){
            scalar_mult_add_wvec( psi+i, p+i, a, psi+i );
            scalar_mult_add_wvec( r+i, mp+i, -a, r+i );
	    rsq += (double)magsq_wvec( r+i );
        }
	g_doublesum( &rsq );
/* if(this_node==0)printf("congrad2: iter %d, rsq %e, pkp %e, a %e\n",
iteration,(double)rsq,(double)pkp,(double)a );  */

	if( rsq <= rsqstop ){
	    *final_rsq_ptr= (Real)rsq;
	    for( i=XUP; i <= TUP; i++) {
		cleanup_gather(tag[i]);
		cleanup_gather(tag[OPP_DIR(i)]);
#ifdef LU
		cleanup_gather(tag2[i]);
		cleanup_gather(tag2[OPP_DIR(i)]);
#endif
	    }
dtime += dclock();
/* if(this_node==0)printf("CONGRAD2: time = %e iters = %d mflops = %e\n",
dtime,iteration,(double)(2840.0*volume*iteration/(1.0e6*dtime*numnodes())) ); */
            cleanup_dslash_temps();
            cleanup_tmp_links();

            FORALLSITES(i,s){
                s->psi = psi[i];
            }
            free(psi); free(chi); free(tmp);
            free(mp); free(p); free(r);

	    return (iteration);
	}

	b = (Real)(rsq/oldrsq);
	FORMYSITESDOMAIN(i,s){
	    scalar_mult_add_wvec( r+i, p+i, b, p+i );
	}

    } while( iteration%niter != 0);

    for( i=XUP; i <= TUP; i++) {
	cleanup_gather(tag[i]);
	cleanup_gather(tag[OPP_DIR(i)]);
#ifdef LU
	cleanup_gather(tag2[i]);
	cleanup_gather(tag2[OPP_DIR(i)]);
#endif
    }

    /*  number of restarts is hard coded */
    if( iteration < CONGRAD_RESTART*niter ) goto start;
    *final_rsq_ptr= (Real)rsq;
    if( rsq > rsqstop ){
	if(this_node==0)printf("No convergence in d_congrad2 size_r= %.2g \n",
                        sqrt(rsq/source_norm));
    }

    cleanup_dslash_temps();
    cleanup_tmp_links();

    FORALLSITES(i,s){
        s->psi = psi[i];
    }
    free(psi); free(chi); free(tmp);
    free(mp); free(p); free(r);

    return(iteration);
}
