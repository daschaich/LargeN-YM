/* Reference:
* Hypercubic smeared links for dynamical fermions.
* By Anna Hasenfratz, Roland Hoffmann, Stefan Schaefer.
* JHEP 0705:029,2007. [hep-lat/0702028]
*/

#include "cl_dyn_includes.h"

void stout_force2(int dir2);
void stout_force3(int dir3, int dir2);
void Sigma_update1 (int dir, su3_matrix_f* sigma_off, su3_matrix_f* stp,
                    su3_matrix_f** l1, Real alpha1, Real alpha2, int sigfresh);
void make_2hermitian_f(su3_matrix_f *A);
void compute_sigma23(su3_matrix_f* sig,
                     su3_matrix_f* lnk1, su3_matrix_f* lnk2,
                     su3_matrix_f* lambda1, su3_matrix_f* lambda2,
                     int dir1, int dir2);
#ifdef SF
void dfatlink_deta(su3_matrix_f* lnk1, su3_matrix_f* lnk2,
                   su3_matrix_f* lambda1, su3_matrix_f* lambda2,
                   int dir1, int dir2);
#endif

/*
* Derivative of the NHYP link wrt to the constinuent links
* Two contributions:
* a) The (thin) link itself : added to Sigma in Sigma_update1
* b) The (fat) staples      : constructed in compute-sigma23 and put in SimgaH
*
*/

/*** first-level force *****************************************************/
void stout_force1()
{

    int dir, dir2;
#if(SMEAR_LEVEL==1)
    int i;
    site* s;
#endif

    /* loop over the link directions, compute sigma from the link itself,
       contruct the Lambdas */

    for (dir = XUP; dir <= TUP; dir++)
    {
  Sigma_update1 (dir, Sigma[dir], Staple3[dir], LambdaU,
                       1.-alpha_smear[0], alpha_smear[0]/6., 0);
    } /*dir */

    /* contruct the field_offsets which point to the links in dir2 direction
       where dir1 is excluded. Here, this makes no sense, however this
       mechanism makes the compute_sigma23 routine re-useable on each level
     */

    for (dir2 = XUP; dir2 <= TUP; dir2++)
    {
  for (dir = XUP; dir <= TUP; dir++)if (dir!=dir2)
  {

#if (SMEAR_LEVEL>1)
      compute_sigma23(SigmaH[dir],
                            hyplink2[dir2][dir], hyplink2[dir][dir2],
                LambdaU[dir], LambdaU[dir2], dir, dir2 );
#ifdef SF
/* in SF, if the original call was from fermion_coupling,            */
/* at this point we compute this level's contribution to K/g^2       */
            if(sf_coupling_flag==SF_COUPLING)
                dfatlink_deta(hyplink2[dir2][dir], hyplink2[dir][dir2],
                  LambdaU[dir], LambdaU[dir2], dir, dir2 );
#endif /* SF */

#else  /* SMEAR_LEVEL==1 */
      compute_sigma23(SigmaH[dir],
                            gauge_field_thin[dir], gauge_field_thin[dir2],
                LambdaU[dir], LambdaU[dir2], dir, dir2 );
            FORALLDYNLINKS(i,s,dir)
            {
               add_su3_matrix_f(Sigma[dir]+i,SigmaH[dir]+i,Sigma[dir]+i);
            }
#ifdef SF
            /* this level's contribution to K/g^2       */
            if(sf_coupling_flag==SF_COUPLING)
                dfatlink_deta(gauge_field_thin[dir], gauge_field_thin[dir2],
                  LambdaU[dir], LambdaU[dir2], dir, dir2 );
#endif /* SF */
#endif /* SMEAR_LEVEL */

  } /* dir */

#if(SMEAR_LEVEL>1)
  stout_force2(dir2);
#endif

    } /* dir2 */

} /* stout_force 1 */


#if (SMEAR_LEVEL>1)
/*** second-level force ****************************************************/
void stout_force2(int dir2) {
  register int i;
  register site *s;
  int dir, dir3, dir4;
#if (SMEAR_LEVEL==3)
  int imap[4][4]={{0,0,1,2},{0,0,0,3},{1,0,0,0},{2,3,0,0}};
  int iimap;
#endif

  /* dir3 is the main direction of the twice smeared hyplink2
     dir2 is the secondary direction of the twice smeared hyplink2
     dir is the main direction of the once-smeared link
     */

  for (dir3 = XUP; dir3 <= TUP; dir3++) if (dir3 != dir2)
  {
    Sigma_update1 (dir3, SigmaH[dir3],Staple2[dir2][dir3],Lambda1,
        (1.-alpha_smear[1]), alpha_smear[1] / 4.,1);
  } /*dir3 */

  for (dir3 = XUP; dir3 <= TUP; dir3++) if( dir3 != dir2)
  {
    for (dir = XUP; dir<= TUP; dir++) if (dir3!=dir && dir!= dir2)
    {
      for (dir4=XUP;dir4<=TUP;dir4++)
      {
        if (dir4!=dir && dir4!=dir2 && dir4 !=dir3) break;
      }

#if (SMEAR_LEVEL==3)
      compute_sigma23(SigmaH[dir],
          hyplink1[dir4][dir], hyplink1[dir4][dir3],
          Lambda1[dir], Lambda1[dir3],dir,dir3 );

#else  /* SMEAR_LEVEL==2 */
      compute_sigma23(SigmaH[dir],
          gauge_field_thin[dir], gauge_field_thin[dir3],
          Lambda1[dir], Lambda1[dir3],dir,dir3 );
      FORALLDYNLINKS(i,s,dir)
      {
        add_su3_matrix_f(Sigma[dir]+i,SigmaH[dir]+i,Sigma[dir]+i);
      }
#endif /* SMEAR_LEVEL */

    } /* dir */

    /* this part is really awkward: The stout_force3 is symmetric in the arguments,
       but the input is not. So one can add the dir2,dir3 and the dir3,dir2 terms.
       to do so, one has to stort them.....I don't like this at all.
       To save some memory, store only the upper triangular part of the 4x4 matrix
       For that, only a 4 'vector' is necessary, which field to use is given by the
       imap array (hard-coded above)
       */

#if (SMEAR_LEVEL==3)
    if (dir2<dir3)
    {
      iimap=imap[dir2][dir3];
      FORALLSITES(i,s)for(dir=XUP;dir<=TUP;dir++)
      {
        su3mat_copy_f(SigmaH[dir]+i,SigmaH2[iimap][dir]+i);
      }

    } else {
      iimap=imap[dir2][dir3];
      FORALLSITES(i,s)for(dir=XUP;dir<=TUP;dir++)
      {
        add_su3_matrix_f(SigmaH[dir]+i,SigmaH2[iimap][dir]+i,
            SigmaH[dir]+i);
      }
      stout_force3(dir2, dir3);
    }
#endif /* SMEAR_LEVEL==3 */

  } /* dir3 */

} /* stout_force 2 */
#endif /* SMEAR_LEVEL>1 */


#if (SMEAR_LEVEL==3)
/*** third-level force *****************************************************/
void stout_force3(int dir3, int dir2)
{
    int dir, dir1,  dir4;
    int i;
    site* s;

    for (dir = XUP; dir <= TUP; dir++) if (dir != dir2 && dir != dir3)
    {
  for (dir4 = XUP; dir4 <= TUP; dir4++)
        {
            if (dir4!=dir && dir4!=dir3 && dir4!=dir2 ) break;
        }
  Sigma_update1 (dir, SigmaH[dir], Staple1[dir4][dir],Lambda2,
                       1.-alpha_smear[2], alpha_smear[2]/2., 1);
    } /*dir */


    for (dir = XUP; dir <= TUP; dir++)  if (dir != dir2 && dir != dir3)
    {
  for (dir1 = XUP; dir1 <= TUP; dir1++)
        {
            if(dir1!=dir && dir1!=dir2 && dir1!=dir3)
      {
          compute_sigma23(tempmat_nhyp1,
                                gauge_field_thin[dir], gauge_field_thin[dir1],
                    Lambda2[dir],Lambda2[dir1],dir,dir1 );
          FORALLDYNLINKS(i,s,dir)
                {
                  add_su3_matrix_f(Sigma[dir]+i,tempmat_nhyp1+i,Sigma[dir]+i);
                }
            }
  }
    }
} /* stout_force 3 */
#endif /* SMEAR_LEVEL==3 */

/*** compute Gamma *********************************************************/
void Sigma_update1 (int dir, su3_matrix_f* sigma_off, su3_matrix_f* stp,
                    su3_matrix_f** lambda, Real alpha1,
                    Real alpha2, int sigfresh)
{
    register int i;
    register site* s;

    su3_matrix_f Gamma, Qisqrt, Q, Omega, tmat, SigmaOmega;
    Real f[NCOL];
    Real bb[NCOL][NCOL];
    complex traces[NCOL];
    complex ctmp;
#if (NCOL>2)
    int j;
    complex ctmp1, ctmp2, ctmp3;
    su3_matrix_f Q2; 
#if (NCOL>3)
    su3_matrix_f Q3,tmat2;
#endif
#endif

    FORALLDYNLINKS(i,s,dir) {
  /* make Omega, Q, Q^2 */
  Omega=stp[i];
  mult_su3_an_f(&Omega,&Omega,&Q);
        /* IR regulator, see clover_xxx/defines.h               */
        /* note that derivatives of Q are unchanged             */
        scalar_add_diag_su3_f(&Q,IR_STAB);

  /* compute inverse sqrt */
#ifndef NHYP_DEBUG
  compute_fhb(&Q, f, bb, 1);
#else
        compute_fhb(&Omega, &Q, f, bb, 1);
#endif

#if (NCOL==2)
        clear_su3mat_f(&Qisqrt);
        Qisqrt.e[0][0].real=f[0];
        Qisqrt.e[1][1].real=f[0];
#else
  mult_su3_nn_f(&Q, &Q, &Q2);
  scalar_mult_su3_matrix_f( &Q, f[1], &tmat);
  scalar_mult_add_su3_matrix_f(&tmat, &Q2, f[2], &Qisqrt);
#if (NCOL>3)
  mult_su3_nn_f(&Q, &Q2, &Q3);
  scalar_mult_add_su3_matrix_f(&Qisqrt, &Q3, f[3], &Qisqrt);
#endif
  scalar_add_diag_su3_f(&Qisqrt, f[0]);
#endif

  /* we'll need Sigma*Omega a few times */
  mult_su3_nn_f(sigma_off+i, &Omega, &SigmaOmega);

  /* now the B matrices and their traces with Sigma*Omega*/
  ctmp=trace_su3_f(&SigmaOmega);
#if (NCOL==2)
  traces[0].real=ctmp.real*bb[0][0];
  traces[0].imag=ctmp.imag*bb[0][0];
  clear_su3mat_f(&Gamma);
  c_scalar_add_diag_su3_f(&Gamma, &traces[0]);
#else
  ctmp1=complextrace_su3_f(&Q,&SigmaOmega);
  ctmp2=complextrace_su3_f(&Q2,&SigmaOmega);
#if (NCOL>3)
  ctmp3=complextrace_su3_f(&Q3,&SigmaOmega);
#endif
  for (j=0;j<NCOL;j++)
  {
      traces[j].real=ctmp.real*bb[j][0]+ctmp1.real*bb[j][1]
                           +ctmp2.real*bb[j][2]
#if (NCOL>3)
         +ctmp3.real*bb[j][3]
#endif
        ;
      traces[j].imag=ctmp.imag*bb[j][0]+ctmp1.imag*bb[j][1]
                           +ctmp2.imag*bb[j][2]
#if (NCOL>3)
         +ctmp3.imag*bb[j][3]
#endif
        ;
  }

  /* the contributions to A tr(B_i Sigma Omega) Q^(i) */
  c_scalar_mult_su3mat_f(&Q, &traces[1], &tmat);
  c_scalar_mult_add_su3mat_f(&tmat , &Q2, &traces[2], &Gamma);
#if (NCOL>3)
  c_scalar_mult_add_su3mat_f(&Gamma, &Q3, &traces[3], &Gamma);
#endif
  c_scalar_add_diag_su3_f(&Gamma, &traces[0]);

  /* the terms propto f_i */
  scalar_mult_add_su3_matrix_f(&Gamma,&SigmaOmega,f[1],&Gamma);
  mult_su3_nn_f(&SigmaOmega,&Q,&tmat);
  scalar_mult_add_su3_matrix_f(&Gamma,&tmat,f[2],&Gamma);
  mult_su3_nn_f(&Q,&SigmaOmega,&tmat);
  scalar_mult_add_su3_matrix_f(&Gamma,&tmat,f[2],&Gamma);
#if (NCOL>3)
  mult_su3_nn_f(&tmat,&Q,&tmat2);
  scalar_mult_add_su3_matrix_f(&Gamma,&tmat2,f[3],&Gamma);
  mult_su3_nn_f(&SigmaOmega,&Q2,&tmat);
  scalar_mult_add_su3_matrix_f(&Gamma,&tmat ,f[3],&Gamma);
  mult_su3_nn_f(&Q2,&SigmaOmega,&tmat);
  scalar_mult_add_su3_matrix_f(&Gamma,&tmat ,f[3],&Gamma);
#endif

#endif
  /* Gamma = (A+A^+)Omega^+ + Q^{-1/2}Sigma */
  make_2hermitian_f(&Gamma);
  mult_su3_na_f(&Gamma,&Omega,&tmat);

  mult_su3_nn_f(&Qisqrt,sigma_off+i,  &Gamma);
  add_su3_matrix_f(&Gamma,&tmat,&Gamma);
  scalar_mult_su3_matrix_f(&Gamma,alpha2,&tmat);
  su3_adjoint_f(&tmat,lambda[dir]+i);

/* the derivative which contributes to the globaln to the new global Sigma
   If this is the first level, then Sigma has to be initiallized. On later
   levels, we accumulate the respecive contributions
*/
        if (sigfresh==0)
            scalar_mult_su3_matrix_f(&Gamma, alpha1, Sigma[dir]+i);
        else
            scalar_mult_add_su3_matrix_f(Sigma[dir]+i, &Gamma, alpha1,
                                   Sigma[dir]+i);
    }/* dynlinks */
} /* Sigma_update1 */


/*** compute next-level Sigma **********************************************/
void compute_sigma23(su3_matrix_f* sig,
                     su3_matrix_f* lnk1, su3_matrix_f* lnk2,
                     su3_matrix_f* lambda1, su3_matrix_f* lambda2,
                     int dir1, int dir2)
{
    register int i;
    register site *st;
    msg_tag *tag0, *tag1, *tag2, *tag3, *tag4;
    su3_matrix_f tmat1, tmat2;

/* construct the next-level Sigma from lnk(s) and lambda(s).
   lambda = Gamma^\dagger
   two staples x three links each = six terms.
   dir1 is the direction of the link and the "main" direction of sigma
   dir2 is the other direction that defines the staple
*/

    /* get link[dir2] from direction dir1   */
    tag0 = start_gather_field( lnk2, sizeof(su3_matrix_f), dir1,
                               EVENANDODD, gen_pt[0] );

    /* get Lambda[dir2] from direction dir1 */
    tag2 = start_gather_field( lambda2, sizeof(su3_matrix_f), dir1,
                               EVENANDODD, gen_pt[2] );

    wait_gather(tag0);
    wait_gather(tag2);

    /* fix boundary values for SF */
    FORALLSITES(i,st) {
        gen_pt[0][i] = CHOOSE_NBR(i,st,dir1,linkf_bndr_up[dir2],0);
        gen_pt[2][i] = CHOOSE_NBR(i,st,dir1,linkf_zero[dir2],2);
    }

    /* get link[dir1] from direction dir2   */
    tag1 = start_gather_field( lnk1, sizeof(su3_matrix_f), dir2,
                               EVENANDODD, gen_pt[1] );

    /* get Lambda[dir1] from direction dir2 */
    tag3 = start_gather_field( lambda1, sizeof(su3_matrix_f), dir2,
                               EVENANDODD, gen_pt[3] );

    /* Lower staple: prepared at x-dir2 and stored in tempmat_nhyp2,
       then gathered to x.
    */

    FORALLSITES(i,st){
        /* "term2" */
        mult_su3_nn_f( lambda1+i, (su3_matrix_f *)gen_pt[0][i], &tmat1 );
        mult_su3_an_f( &tmat1, lnk2+i, tempmat_nhyp2+i );
        /* "term3" */
        mult_su3_nn_f( lnk1+i, (su3_matrix_f *)gen_pt[2][i], &tmat1 );
        mult_su3_an_f( &tmat1, lnk2+i, &tmat2 );
        add_su3_matrix_f(tempmat_nhyp2+i, &tmat2, tempmat_nhyp2+i );
        /* "term4" */
        mult_su3_nn_f( lnk1+i, (su3_matrix_f *)gen_pt[0][i], &tmat1 );
        mult_su3_an_f( &tmat1, lambda2+i, &tmat2 );
        add_su3_matrix_f(tempmat_nhyp2+i, &tmat2, tempmat_nhyp2+i );
    }

    /* gather staple from direction -dir2 to "home" site */
    tag4 = start_gather_field( tempmat_nhyp2, sizeof(su3_matrix_f),
                               OPP_DIR(dir2), EVENANDODD, gen_pt[4] );

    wait_gather(tag1);
    wait_gather(tag3);

    /* fix boundary values for SF */
    FORALLSITES(i,st) {
        gen_pt[1][i] = CHOOSE_NBR(i,st,dir2,linkf_bndr_up[dir1],1);
        gen_pt[3][i] = CHOOSE_NBR(i,st,dir2,linkf_zero[dir1],3);
    }


    /* Upper staple */
    FORALLDYNLINKS(i,st,dir1){
        /* "term1" */
        mult_su3_nn_f( lambda2+i, (su3_matrix_f *)gen_pt[1][i], &tmat1 );
        mult_su3_na_f( (su3_matrix_f *)gen_pt[0][i], &tmat1, sig+i );
        /* "term5" */
        mult_su3_na_f( (su3_matrix_f *)gen_pt[2][i],
                       (su3_matrix_f *)gen_pt[1][i], &tmat1 );
        mult_su3_na_f( &tmat1, lnk2+i, &tmat2 );
        add_su3_matrix_f(sig+i, &tmat2, sig+i);
        /* "term6" */
        mult_su3_na_f( (su3_matrix_f *)gen_pt[0][i],
                       (su3_matrix_f *)gen_pt[3][i], &tmat1 );
        mult_su3_na_f( &tmat1, lnk2+i, &tmat2 );
        add_su3_matrix_f(sig+i, &tmat2, sig+i);
    }

    /* finally add the lower staple. */
    wait_gather(tag4);

    FORALLDYNLINKS(i,st,dir1){
  add_su3_matrix_f( sig+i, (su3_matrix_f *)gen_pt[4][i], sig+i );
    }

    cleanup_gather(tag0);
    cleanup_gather(tag1);
    cleanup_gather(tag2);
    cleanup_gather(tag3);
    cleanup_gather(tag4);

} /* compute_sigma23 */




/*** helper ****************************************************************/
void make_2hermitian_f(su3_matrix_f *A) {
    int i,j;
    for (i=0;i<NCOL;i++) for (j=i;j<NCOL;j++)
    {
  A->e[i][j].real=(A->e[i][j].real+A->e[j][i].real);
  A->e[i][j].imag=(A->e[i][j].imag-A->e[j][i].imag);
  A->e[j][i].real= A->e[i][j].real;
  A->e[j][i].imag=-A->e[i][j].imag;
    }
}
