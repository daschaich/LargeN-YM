/********   sf_make_boundary.c   ***************************************
************************************************************************
modified from milc-7.4.0 generic_schroed/make_schroed_lattice.c
************************************************************************

Handling Schroedinger Functional (SF) boundary conditions:

The gather mechanism assumes periodic boundary conditions,
and we will have to supersede it if required by SF.

Boundary values structure (from lattice.h):

EXTERN  su3_matrix   link_bndr_up[3],   link_driv_up[3],
                     link_bndr_dn[3],   link_driv_dn[3];
EXTERN  su3_matrix_f linkf_bndr_up[3],  linkf_driv_up[3],
                     linkf_bndr_dn[3],  linkf_driv_dn[3];

The t=0 boundary values for link and linkf are also placed into *lattice .
Macro CHOOSE_NBR supersedes the gather mechanism
when we need to gather a link from the boundary at "t=nt" .
Because we do not update the spatial links at t=0,
we don't care what would-be force is computed for them.

The dynamical fermions live on 1 <= t <= nt-1 .
Hence, at t=0, all wilson_vectors and [functions of] the clover term
must be zero prior to any gather.  Typically (e.g. in dslash) we gather
from all sites, but use the gathered values only on t>0 sites.

Employ "twisted" boundary conditions in spatial directions for the fermions.

************************************************************************
*** set_boundary_values() **********************************************

  set up boundary values of gauge field:
  global linkf and link, for bndr and for driv
  local t=0 boundary links also placed into *lattice
  all wilson vectors set to zero at t=0.
  Also initialize ferm_phases.

  bc_flag  --  choice between coded boundary values, current options:

  0  TRIVIAL    U_k=1, for debugging only
  1  ALPHA      SU(2) only, alpha collaboration's favorite
  2  ALPHA_ONE  SU(3) only, used by alpha collab.
  3  ALPHA_TWO  SU(3) only, also used by alpha collab.
  4  SF_SYM     SU(4) only for now; analogue of ALPHA_ONE in that
                boundary values situated symmetrically in fund. domain

  c_t -- multiplies boundary linkf but not boundary link.
  Implements coefficient of space-time boundary plaquettes in gauge action.
  ct_flag's current options:

  0  TREE_LEVEL   c_t = 1
  1  ONE_LOOP     c_t = 1 + c^(1) g_0^2 (see below for numbers)

BACKGROUND_FIELD_X -- flag in defines.h

  0 turns off background field in x direction; requires ny=nz
  1 turns on background field in x direction; requires nx=ny=nz

ferm_twist_phase is always applied in all three spatial directions.

***********************************************************************/


#include "cl_dyn_includes.h"


void set_boundary_values(){
register int j,k,dir;
Real factor[NCOL],b0r[NCOL],b0i[NCOL],btr[NCOL],bti[NCOL];
Real angle;

/* for comparison to alpha collab *********************************
   the following supersedes the input values

Real gsq=6./beta;

  clov_c = (1+gsq*(-0.454+gsq*(-0.175+gsq*(0.012+gsq*0.045))))
           / (1-gsq*0.72);

  kappa = 0.125
    + gsq*(0.008439857+gsq*(0.0085+gsq*(-0.0272+gsq*(0.0420-gsq*0.0204))));

  u0=1.;

  if(this_node==0)printf("PARAMETERS clov_c= %f  kappa= %f  u0= %f  ",
                         clov_c, kappa, u0 );
*/

  if(this_node==0)printf("SF PARAMETERS ");

/*****************************************************************/

   switch(ct_flag){
/*  tree-level value  */
   case TREE_LEVEL:
      c_t = 1.0;
      break;
#if(NCOL==3)
#if(DIMF==3)
/*  one-loop value    */
   case ONE_LOOP:
      c_t = 1.0 - 0.3042/beta;  /* with two fundamental quarks */
/*    c_t = 1.0 - 0.534/beta;     quencehd */
      break;
/* two-loop value
   two-loop correction:  -0.029 + 2*0.002 = -0.025
*/
   case TWO_LOOP:
      c_t = 1.0 - 0.3042/beta - 0.025*(6./beta)*(6./beta) ;
      break;
#endif
#endif
   default:
      if(this_node==0)printf("Illegal value of ct_flag, internal code %d\n", ct_flag);
      terminate(1);
   }
   if(this_node==0)printf("c_t= %f\n", c_t );

/* Check that geometry is consistent with value of BACKGROUND_FIELD_X */

   if(BACKGROUND_FIELD_X==0){
      if(ny!=nz){
         if(this_node==0)printf("ny=nz required!\n");
         terminate(1);
      }
   }

   if(BACKGROUND_FIELD_X==1){
      if(nx!=ny || nx!=nz){
         if(this_node==0)printf("nx=ny=nz required!\n");
         terminate(1);
      }
   }

/* prepare fermion phases */

   ferm_phases[0].real = cos( (double)(ferm_twist_phase*PI/(Real)nx) );
   ferm_phases[0].imag = sin( (double)(ferm_twist_phase*PI/(Real)nx) );
   ferm_phases[1].real = cos( (double)(ferm_twist_phase*PI/(Real)ny) );
   ferm_phases[1].imag = sin( (double)(ferm_twist_phase*PI/(Real)ny) );
   ferm_phases[2].real = cos( (double)(ferm_twist_phase*PI/(Real)nz) );
   ferm_phases[2].imag = sin( (double)(ferm_twist_phase*PI/(Real)nz) );

/* set boundary values of gauge field */

/* create boundary values. */

   switch(bc_flag){
/* trivial boundary values (for debugging)        */
      case TRIVIAL:
	 for(j=0; j<NCOL; j++){
            b0r[j] = 1.;
            b0i[j] = 0.;
            btr[j] = 1.;
            bti[j] = 0.;
	 }
	 break;
/* alpha collaboration's favorite choices of boundary values */
#if (NCOL==2)
      case ALPHA:
	    angle = PI/(Real)(4*nt);
            b0r[0] = cos((double)angle);
            b0i[0] = sin((double)angle);
            b0r[1] = b0r[0];
            b0i[1] = -b0i[0];
            btr[0] = cos((double)(3*angle));
            bti[0] = sin((double)(3*angle));
            btr[1] = btr[0];
            bti[1] = -bti[0];
	    break;
#elif (NCOL==3)
      case ALPHA_ONE:
	    angle = PI/(Real)(3*nt);
            b0r[2] = cos((double)angle);
            b0i[2] = sin((double)angle);
	    b0r[1] = 1.;
	    b0i[1] = 0.;
	    b0r[0] = b0r[2];
            b0i[0] = -b0i[2];
            btr[2] = cos((double)(2.0*angle));
            bti[2] = sin((double)(2.0*angle));
            btr[1] = b0r[2];
            bti[1] = b0i[2];
            btr[0] = cos((double)(3.0*angle));
            bti[0] = -sin((double)(3.0*angle));
	    break;
      case ALPHA_TWO:
	    angle = PI/(Real)(6*nt);
            b0r[2] = cos((double)angle);
            b0i[2] = sin((double)angle);
	    b0r[1] = 1.;
	    b0i[1] = 0.;
	    b0r[0] = b0r[2];
            b0i[0] = -b0i[2];
            btr[2] = cos((double)(3.0*angle));
            bti[2] = sin((double)(3.0*angle));
            btr[1] = cos((double)(2.0*angle));
            bti[1] = sin((double)(2.0*angle));
            btr[0] = cos((double)(5.0*angle));
            bti[0] = -sin((double)(5.0*angle));
	    break;
#elif (NCOL==4)
/*
Biagio Lucini, Gregory Moraitis, Phys.Lett.B668:226-232,2008.
e-Print: arXiv:0805.2913
*/
      case SF_SYM:
	    angle = sqrt(2)*PI/(Real)(4*nt);
	    b0r[3] = cos((double)angle);
	    b0i[3] = sin((double)angle);
            b0r[2] = cos((double)(PI/(Real)(2*nt) - angle));
            b0i[2] = sin((double)(PI/(Real)(2*nt) - angle));
	    b0r[1] = b0r[2];
	    b0i[1] = -b0i[2];
	    b0r[0] = b0r[3];
            b0i[0] = -b0i[3];
	    btr[3] = cos((double)(PI/(Real)(2*nt) + angle));
	    bti[3] = sin((double)(PI/(Real)(2*nt) + angle));
            btr[2] = cos((double)(PI/(Real)nt - angle));
            bti[2] = sin((double)(PI/(Real)nt - angle));
            btr[1] = btr[2];
            bti[1] = -bti[2];
            btr[0] = btr[3];
            bti[0] = -bti[3];
	    break;
#else
        if(this_node==0)printf("No boundary values for NCOL>4 !\n");
        terminate(1);
#endif
      default:
	 if(this_node==0)printf("Bad bc_flag, internal code %d\n",bc_flag);
	 terminate(1);
   }

/* factors for derivatives (all dirs)
   linkf_driv_dn.e[k][k] = i*factor[k]*linkf_bndr_dn.e[k][k]
   linkf_driv_up.e[k][k] = -i*factor[k]*linkf_bndr_up.e[k][k]
*/

#if (NCOL==2)
   factor[0] = 1.0/(Real)nt;
   factor[1] = -factor[0];
#elif (NCOL==3)
   factor[0] = 1.0/(Real)nt;
   factor[1] = -0.5*factor[0];
   factor[2] = factor[1];
#else /* SU(4) */
   factor[0] = -0.5/(Real)nt;
   factor[1] = -0.5/(Real)nt;
   factor[2] = 0.5/(Real)nt;
   factor[3] = 0.5/(Real)nt;
#endif

/* the rest works for any NCOL      */
/* boundary values for linkf        */

   for(dir=XUP;dir<=ZUP;dir++){
#ifdef NHYP
      clear_su3mat_f(linkf_zero+dir);
#endif

/*  Check if we want background field in x direction  */
      if(dir==XUP && BACKGROUND_FIELD_X==0){
         for(j=0; j<NCOL; j++) for(k=0; k<NCOL; k++){
            if(j!=k){
               linkf_bndr_dn[dir].e[j][k] = cmplx(0.0,0.0);
               linkf_bndr_up[dir].e[j][k] = cmplx(0.0,0.0);
               linkf_driv_dn[dir].e[j][k] = cmplx(0.0,0.0);
               linkf_driv_up[dir].e[j][k] = cmplx(0.0,0.0);
            }
            else{
               linkf_bndr_dn[dir].e[j][k].real = 1.;
               linkf_bndr_dn[dir].e[j][k].imag = 0.;
               linkf_bndr_up[dir].e[j][k].real = 1.;
               linkf_bndr_up[dir].e[j][k].imag = 0.;
               linkf_driv_dn[dir].e[j][k].real = 0.;
               linkf_driv_dn[dir].e[j][k].imag = 0.;
               linkf_driv_up[dir].e[j][k].real = 0.;
               linkf_driv_up[dir].e[j][k].imag = 0.;
            }
         }
      }
      else{
         for(j=0; j<NCOL; j++) for(k=0; k<NCOL; k++){
            if(j!=k){
               linkf_bndr_dn[dir].e[j][k] = cmplx(0.0,0.0);
               linkf_bndr_up[dir].e[j][k] = cmplx(0.0,0.0);
               linkf_driv_dn[dir].e[j][k] = cmplx(0.0,0.0);
               linkf_driv_up[dir].e[j][k] = cmplx(0.0,0.0);
            }
            else{
               linkf_bndr_dn[dir].e[j][k].real = b0r[j];
               linkf_bndr_dn[dir].e[j][k].imag = b0i[j];
               linkf_bndr_up[dir].e[j][k].real = btr[j];
               linkf_bndr_up[dir].e[j][k].imag = bti[j];
               linkf_driv_dn[dir].e[j][k].real = -factor[j]*b0i[j];
               linkf_driv_dn[dir].e[j][k].imag = factor[j]*b0r[j];
               linkf_driv_up[dir].e[j][k].real = factor[j]*bti[j];
               linkf_driv_up[dir].e[j][k].imag = -factor[j]*btr[j];
            }
         }
      }

/* boundary values in the fermion's irrep                         */
/* boundary link from linkf, then twist phase                     */
      make_fermion_rep_matrix( linkf_bndr_dn+dir, link_bndr_dn+dir );
      sf_phase( link_bndr_dn+dir, dir );
      make_fermion_rep_matrix( linkf_bndr_up+dir, link_bndr_up+dir );
      sf_phase( link_bndr_up+dir, dir );
/* the parametric derivatives, according to the Leibniz rule,
   inputs are linkf_bndr_xx and factor (not linkf_driv_xx !!)     */
      make_fermion_rep_driv( linkf_bndr_dn+dir, link_driv_dn+dir,
                             factor, PLUS );
      sf_phase( link_driv_dn+dir, dir );
      make_fermion_rep_driv( linkf_bndr_up+dir, link_driv_up+dir,
                             factor, MINUS );
      sf_phase( link_driv_up+dir, dir );
   }

/* now that link_bndr etc are set, multiply linkf_bndr etc by c_t */

   for(dir=XUP;dir<=ZUP;dir++){
      for(j=0; j<NCOL; j++) {
         linkf_bndr_dn[dir].e[j][j].real *= c_t;
         linkf_bndr_dn[dir].e[j][j].imag *= c_t;
         linkf_bndr_up[dir].e[j][j].real *= c_t;
         linkf_bndr_up[dir].e[j][j].imag *= c_t;
         linkf_driv_dn[dir].e[j][j].real *= c_t;
         linkf_driv_dn[dir].e[j][j].imag *= c_t;
         linkf_driv_up[dir].e[j][j].real *= c_t;
         linkf_driv_up[dir].e[j][j].imag *= c_t;
      }
   }
}  /* set_boundary_values */


/***********************************************************************/
void initialize_sf_lattice(){
register int i,j,k,dir;
register site *s;
#ifdef NHYP
register int otherdir;
#endif

   FORALLSITES(i,s){
      if(s->t == 0){
/* set copies of boundary links in *lattice  */
	 for(dir=XUP;dir<=ZUP;dir++){
            s->linkf[dir]=linkf_bndr_dn[dir];
            s->link[dir]=link_bndr_dn[dir];
#ifdef NHYP
/* likewise we copy the spatial boundary values into the fields
   containing the fat and thin link, the Omega's and the V's.
   For tensors used in the force calculation we set
   the spatial boundary values to zero.
   SigmaH2 requires no initialization.
*/
            gauge_field_thin[dir][i]=linkf_bndr_dn[dir];
            gauge_field[dir][i]=linkf_bndr_dn[dir];
            Staple3[dir][i]=linkf_bndr_dn[dir];
            clear_su3mat_f(Sigma[dir]+i);
            clear_su3mat_f(SigmaH[dir]+i);
            clear_su3mat_f(LambdaU[dir]+i);
#if (SMEAR_LEVEL>1)
            clear_su3mat_f(Lambda1[dir]+i);
#endif
#if (SMEAR_LEVEL==3)
            clear_su3mat_f(Lambda2[dir]+i);
#endif
            for(otherdir=XUP;otherdir<=TUP;otherdir++){
               if(otherdir!=dir){
#if (SMEAR_LEVEL==3)
                  Staple1[otherdir][dir][i]=linkf_bndr_dn[dir];
                  hyplink1[otherdir][dir][i]=linkf_bndr_dn[dir];
#endif
#if (SMEAR_LEVEL>1)
                  Staple2[otherdir][dir][i]=linkf_bndr_dn[dir];
                  hyplink2[otherdir][dir][i]=linkf_bndr_dn[dir];
#endif
               }
            }
#endif
	 }
/* set all wilson vectors to zero on t=0=nt */
         for(k=0;k<4;k++)for(j=0;j<DIMF;j++){
            s->g_rand.d[k].c[j] = cmplx(0.0,0.0);
            s->psi[0].d[k].c[j]    = cmplx(0.0,0.0);
            s->psi[1].d[k].c[j]    = cmplx(0.0,0.0);
            s->chi[0].d[k].c[j]    = cmplx(0.0,0.0);
            s->chi[1].d[k].c[j]    = cmplx(0.0,0.0);
	    s->tmp.d[k].c[j]    = cmplx(0.0,0.0);
	    s->mp.d[k].c[j]     = cmplx(0.0,0.0);
	    s->p.d[k].c[j]      = cmplx(0.0,0.0);
	    s->r.d[k].c[j]      = cmplx(0.0,0.0);
         }
      }
   }

//   if(this_node==0)printf("SF lattice initialized\n");
}  /* initialize_sf_lattice */
