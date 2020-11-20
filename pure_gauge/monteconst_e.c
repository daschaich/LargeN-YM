// -----------------------------------------------------------------
// Kennedy--Pendleton quasi-heat bath (qhb) on SU(2) subgroups
#include "pg_includes.h"
#define INC 1.0e-10

double action1();

void monteconst_e(int NumStp, double Eint, double delta, double a) {
  register int dir, i;
  register site *st;
  int NumTrj, Nhit, subgrp, ina, inb, ii, parity, count;
  int j, k, kp, cr, nacd, test, index_a[N_OFFDIAG], index_b[N_OFFDIAG];
  Real xr1, xr2, xr3, xr4;
  Real a0 = 0, a1, a2, a3;
  Real v0, v1, v2, v3, vsq;
  Real h0, h1, h2, h3;
  Real r, r2, rho, z;
  Real al, d, xl, xd;
  Real pi2, b3;
  su2_matrix h;
  matrix_f action;
  matrix_f oldlinkvalue;
  matrix_f oldvalue;
  //matrix_f oldvalueneg;
  matrix_f newvalue;
  matrix_f tracematrix;
  
  complex trace;
  double realtrace;
  double Energydiff;
  double Energy;
  double Energyref;
  
  
  
  //double betareference = 9.6;

  Nhit = (int)N_OFFDIAG;    // NCOL * (NCOL - 1) / 2
  pi2 = 2.0 * PI;
  b3 = beta * a* one_ov_N;

  // Set up SU(2) subgroup indices [a][b], always with a < b
  count = 0;
  for (i = 0; i < NCOL; i++) {
    for (j = i + 1; j < NCOL; j++) {
      index_a[count] = i;
      index_b[count] = j;
      count++;
    }
  }
  if (count != Nhit) {      // Sanity check
    node0_printf("ERROR: %d rather than %d subgroups found", count, Nhit);
    terminate(1);
  }
  Energy = action1();
  Energyref = action1();
  /* fix bug by adding loop over NumTrj; before 1 (and only 1) heat bath
     hit was done, regardless of NumStp    */
  for (NumTrj = 0; NumTrj < NumStp; NumTrj++) {
    /* fix bug by looping over odd AND even parity */
    for (parity=ODD;parity<=EVEN;parity++) {
      FORALLUPDIR(dir) {
        // Compute the gauge force
        dsdu_qhb(dir, parity);

        // Now for the qhb updating, looping over SU(2) subgroups
        for (subgrp = 0; subgrp < Nhit; subgrp++) {
          kp=0;
          cr=0;

          // Pick out this SU(2) subgroup
          ina = index_a[subgrp];
          inb = index_b[subgrp];
          FORSOMEPARITY(i, st, parity) {
            // st = &(lattice.[i])
            //scalar_mult_mat(&(st->linkf[dir]), 1, &oldlinkvalue);
            mat_copy_f(&(lattice[i].linkf[dir]), &oldlinkvalue);
            mult_na_f(&(st->linkf[dir]), &(st->staple), &action);
            //mult_na_f(&(st->linkf[dir]), &(st->staple), &oldvalue);
            mult_an_f(&(lattice[i].linkf[dir]), &(lattice[i].staple), &oldvalue);
            //scalar_mult_mat_f(&oldvalue, -1, &oldvalueneg);
            
            /*decompose the action into SU(2) subgroups using Pauli matrix expansion */
            /* The SU(2) hit matrix is represented as v0 + i * Sum j (sigma j * vj)*/
            v0 = action.e[ina][ina].real + action.e[inb][inb].real;
            v3 = action.e[ina][ina].imag - action.e[inb][inb].imag;
            v1 = action.e[ina][inb].imag + action.e[inb][ina].imag;
            v2 = action.e[ina][inb].real - action.e[inb][ina].real;

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
            if ((1.00 - 0.5 * d) > xr4 * xr4)
              nacd=1;

            // Kennedy--Pendleton algorithm
            if (nacd == 0 && al > 2.0) {
              test=0;
              for (k=0;k<20 && test == 0;k++) {
                kp++;
                /*  get four random numbers (add a small increment to prevent taking log(0.)*/
                xr1 = myrand(&(st->site_prn));
                xr1 = log((double)(xr1 + INC));

                xr2 = myrand(&(st->site_prn));
                xr2 = log((double)(xr2 + INC));

                xr3 = myrand(&(st->site_prn));
                xr3 = cos((double)pi2 * xr3);
                d = -(xr2 + xr1 * xr3 * xr3) / al;

                xr4 = myrand(&(st->site_prn));
                if ((1.0 - 0.5 * d) > xr4 * xr4)
                  test = 1;
              }
              if (test != 1)
                node0_printf("site took 20 kp hits\n");
            }

            if (nacd == 0 && al <= 2.0) {/* creutz algorithm */
              cr++;
              xl = exp((double)(-2.0 * al));
              xd = 1.0 - xl;
              test = 0;
              for (k = 0; k < 20 && test == 0; k++) {
                // Get two random numbers
                xr1 = myrand(&(st->site_prn));
                xr2 = myrand(&(st->site_prn));

                r = xl + xd * xr1;
                a0 = 1.00 + log((double)r) / al;
                if ((1.0 - a0 * a0) > xr2 * xr2)
                  test = 1;
              }
              d = 1.0 - a0;
              if (test != 1)
                node0_printf("site took 20 Creutz hits\n");
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
            rho = sqrt((double)rho);

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
            
            mult_an_f(&(st->linkf[dir]), &(st->staple), &newvalue);
            sub_mat_f(&newvalue,&oldvalue,&tracematrix);
            (trace).real=0;
            (trace).imag=0;
            trace_sum_f(&tracematrix, &trace);
            realtrace=(trace).real;
            
            Energydiff = -beta*realtrace;
            Energy = Energy + Energydiff/3.0;
            //Energyref = action1();
            //node0_printf("energy %.8g %.8g %.8g\n",
            //     Energyref, Energy, Energydiff/3.0);
            
            if((Energy>=Eint && Energy <= (Eint+delta))==false)
            {
              mat_copy_f(&oldlinkvalue,&(lattice[i].linkf[dir]));
              Energy = Energy - Energydiff/3.0;
            }    
            
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
        }
      }
    }
   // node0_printf("energy %.8g %.8g %.8g\n",
   //       Energyref, Energy, Energydiff/3.0);
  }
}
double action1() {
  double ssplaq, stplaq, g_act;
  //double betareference=9.6;
  plaquette(&ssplaq, &stplaq);
  g_act = -beta * volume * (ssplaq + stplaq);
  
  return g_act;
}

// -----------------------------------------------------------------
