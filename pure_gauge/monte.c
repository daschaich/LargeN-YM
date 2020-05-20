// -----------------------------------------------------------------
// Kennedy--Pendleton quasi-heat bath (qhb) on SU(2) subgroups
#include "pg_includes.h"
#define INC 1.0e-10

void monte(int NumStp) {
  register int dir,i;
  register site *st;
  int NumTrj, Nhit, index1, ina, inb, ii, parity;
  Real xr1, xr2, xr3, xr4;
  Real a0 = 0, a1, a2, a3;
  Real v0, v1, v2, v3, vsq;
  Real h0, h1, h2, h3;
  Real r, r2, rho, z;
  Real al, d, xl, xd;
  int  k, kp, cr, nacd, test;
  Real pi2, b3;
  su2_matrix h;
  matrix_f action;

  Nhit = 3;
  pi2 = 2.0 * PI;
  b3 = beta / 3.0; // TODO: NCOL or dim*(dim-1)/2?

  /* fix bug by adding loop over NumTrj; before 1 (and only 1) heat bath
     hit was dome, regardless of NumStp    */
  for (NumTrj = 0 ; NumTrj < NumStp; NumTrj++) {
    /* fix bug by looping over odd AND even parity */
    for (parity=ODD;parity<=EVEN;parity++) {
      FORALLUPDIR(dir) {
        /* compute the gauge force */
        dsdu_qhb(dir,parity);
        /* now for the qhb updating */
        for (index1=0;index1<Nhit;index1++) {
          kp=0;
          cr=0;

          /*  pick out an SU(2) subgroup */
          ina=(index1+1) % NCOL;
          inb=(index1+2) % NCOL;
          if (ina > inb) { ii=ina; ina=inb; inb=ii;}

          FORSOMEPARITY(i, st, parity) {
            mult_na_f(&(st->linkf[dir]), &(st->staple), &action);

            /*decompose the action into SU(2) subgroups using Pauli matrix expansion */
            /* The SU(2) hit matrix is represented as a0 + i * Sum j (sigma j * aj)*/
            v0 =  action.e[ina][ina].real + action.e[inb][inb].real;
            v3 =  action.e[ina][ina].imag - action.e[inb][inb].imag;
            v1 =  action.e[ina][inb].imag + action.e[inb][ina].imag;
            v2 =  action.e[ina][inb].real - action.e[inb][ina].real;

            vsq = v0*v0 + v1*v1 + v2*v2 + v3*v3;
            z = sqrt((double)vsq);
            /* Normalize   u */
            v0 = v0/z; v1 = v1/z; v2 = v2/z; v3 = v3/z;

            /* end norm check--trial SU(2) matrix is a0 + i a(j)sigma(j)*/

            /* test 
               if (this_node == 0)printf("v= %e %e %e %e\n",v0,v1,v2,v3);
               if (this_node == 0)printf("z= %e\n",z);
               */
            /* now begin qhb */
            /* get four random numbers */

            /*  get four random numbers (add a small increment to prevent taking log(0.)*/
            xr1 = myrand(&(st->site_prn));
            xr1 = (log((double)(xr1 + INC)));

            xr2 = myrand(&(st->site_prn));
            xr2 = (log((double)(xr2 + INC)));

            xr3 = myrand(&(st->site_prn));
            xr4 = myrand(&(st->site_prn));

            xr3 = cos((double)pi2 * xr3);
            /*
               node0_printf("rand= %e %e %e %e\n", xr1, xr2, xr3, xr4);
               */

            /*
               generate a0 component of su3 matrix

               first consider generating an su(2) matrix h
               according to exp(bg/3 * re tr(h*s))
               rewrite re tr(h*s) as re tr(h*v)z where v is
               an su(2) matrix and z is a real normalization constant 
               let v = z*v. (z is 2*xi in Kennedy--Pendleton notation)
               v is represented in the form v(0) + i*sig*v (sig are pauli)
               v(0) and vector v are real

               let a = h*v and now generate a
               rewrite beta/3 * re tr(h*v) * z as al*a0
               a0 has prob(a0) = n0 * sqrt(1 - a0**2) * exp(al * a0)
               */
            al=b3*z;
            /*if (this_node == 0)printf("al= %e\n",al);*/

            /*
               let a0 = 1 - del**2
               get d = del**2
               such that prob2(del) = n1 * del**2 * exp(-al*del**2)
               */

            d= -(xr2  + xr1*xr3*xr3)/al;

            /*     monte carlo prob1(del) = n2 * sqrt(1 - 0.5*del**2)
                   then prob(a0) = n3 * prob1(a0)*prob2(a0)
                   */

            /* now  beat each  site into submission */
            nacd = 0;
            if ((1.00 - 0.5*d) > xr4*xr4)
              nacd=1;

            // Kennedy--Pendleton algorithm
            if (nacd == 0 && al > 2.0) {
              test=0;
              for (k=0;k<20 && test == 0;k++) {
                kp++;
                /*  get four random numbers (add a small increment to prevent taking log(0.)*/
                xr1=myrand(&(st->site_prn));
                xr1 =  (log((double)(xr1+ INC)));

                xr2=myrand(&(st->site_prn));
                xr2 =  (log((double)(xr2+ INC)));

                xr3=myrand(&(st->site_prn));
                xr4=myrand(&(st->site_prn));

                xr3=cos((double)pi2*xr3);

                d = -(xr2 + xr1*xr3*xr3)/al;
                if ((1.00 - 0.5*d) > xr4*xr4)
                  test = 1;
              }
              if (this_node == 0 && test !=1)
                printf("site took 20 kp hits\n");
            } /* endif nacd */

            if (nacd == 0 && al <= 2.0) {/* creutz algorithm */
              cr++;
              xl = exp((double)(-2.0 * al));
              xd = 1.0 - xl;
              test = 0;
              for (k = 0; k < 20 && test == 0; k++) {
                /*        get two random numbers */
                xr1 = myrand(&(st->site_prn));
                xr2 = myrand(&(st->site_prn));

                r = xl + xd*xr1; 
                a0 = 1.00 + log((double)r) / al;
                if ((1.0 -a0*a0) > xr2*xr2)
                  test = 1;
              }
              d = 1.0 - a0;
              if (this_node == 0 && test !=1) 
                printf("site took 20 creutz hits\n");
            } /* endif nacd */

            /*  generate full su(2) matrix and update link matrix*/

            /* find a0  = 1 - d*/
            a0 = 1.0 - d;
            /* compute r */
            r2 = 1.0 - a0*a0;
            r2 = fabs((double)r2);
            r = sqrt((double)r2);

            /* compute a3 */
            a3=(2.0*myrand(&(st->site_prn)) - 1.0)*r;

            /* compute a1 and a2 */
            rho = r2 - a3*a3;
            rho = fabs((double)rho);
            rho= sqrt((double)rho);

            /*xr2 is a random number between 0 and 2*pi */
            xr2=pi2*myrand(&(st->site_prn));

            a1= rho*cos((double)xr2);
            a2= rho*sin((double)xr2);

            /* now do the updating.  h = a*v^dagger, new u = h*u */
            h0 = a0*v0 + a1*v1 + a2*v2 + a3*v3;
            h1 = a1*v0 - a0*v1 + a2*v3 - a3*v2;
            h2 = a2*v0 - a0*v2 + a3*v1 - a1*v3;
            h3 = a3*v0 - a0*v3 + a1*v2 - a2*v1;

            /* Elements of SU(2) matrix */
            h.e[0][0] = cmplx( h0, h3);
            h.e[0][1] = cmplx( h2, h1);
            h.e[1][0] = cmplx(-h2, h1);
            h.e[1][1] = cmplx( h0,-h3);

            /* update the link */
            left_su2_hit_n_f(&h, ina, inb, &(st->linkf[dir]));
          }
          /* diagnostics
             {Real avekp, avecr;
             avekp=(Real)kp / (Real)(nx*ny*nz*nt/2);
             avecr=(Real)cr / (Real)(nx*ny*nz*nt/2);
             if (this_node ==0)
             printf(" ave kp steps = %e, ave creutz steps = %e\n",
             (double)avekp,(double)avecr);
             }
             */
        } /*  hits */
      } /*  direction */
    }} /* parity and NumTrj */
}
// -----------------------------------------------------------------
