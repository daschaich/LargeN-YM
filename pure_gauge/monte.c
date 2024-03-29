// -----------------------------------------------------------------
// Kennedy--Pendleton quasi-heatbath (qhb) on SU(2) subgroups
#include "pg_includes.h"
#define INC 1.0e-10
//#define DEBUG_PRINT

void monte() {
  register int dir, i;
  register site *s;
  int istep, Nhit, subgrp, ina, inb, parity, count;
  int k, kp, cr, nacd, test, index_a[N_OFFDIAG], index_b[N_OFFDIAG];
  Real xr1, xr2, xr3, xr4;
  Real a0 = 0, a1, a2, a3;
  Real v0, v1, v2, v3, vsq;
  Real h0, h1, h2, h3;
  Real r, r2, rho, z, norm;
  Real al, d, xl, xd, b3 = beta * one_ov_N;
  su2_matrix h;
  matrix action;
#ifdef DEBUG_PRINT
  double ss_plaq, st_plaq, check;
#endif

  // Set up SU(2) subgroup indices [a][b] with a < b
  count = 0;
  Nhit = (int)N_OFFDIAG;    // NCOL * (NCOL - 1) / 2
  for (ina = 0; ina < NCOL - 1; ina++) {
    for (inb = ina + 1; inb < NCOL; inb++) {
      index_a[count] = ina;
      index_b[count] = inb;
      count++;
    }
  }
  if (count != Nhit) {      // Sanity check
    node0_printf("ERROR: %d rather than %d subgroups found", count, Nhit);
    terminate(1);
  }

  // Loop over quasi-heatbath sweeps
  for (istep = 0; istep < qhb_steps; istep++) {
    for (parity = ODD; parity <= EVEN; parity++) {
      FORALLUPDIR(dir) {
        // Compute the gauge force (updating s->staple)
        dsdu_qhb(dir, parity);

        // Now for the qhb updating, looping over SU(2) subgroups
        for (subgrp = 0; subgrp < Nhit; subgrp++) {
          kp = 0;
          cr = 0;

          // Pick out this SU(2) subgroup
          ina = index_a[subgrp];
          inb = index_b[subgrp];
          FORSOMEPARITY(i, s, parity) {
            // Decompose the action into SU(2) subgroups
            // using Pauli matrix expansion
            // The SU(2) hit matrix is represented as
            //   v0 + i * Sum j (sigma j * vj)
            mult_na(&(s->link[dir]), &(s->staple), &action);
#ifdef DEBUG_PRINT
            v0 = action.e[ina][ina].real + action.e[inb][inb].real;
            v3 = action.e[ina][ina].imag - action.e[inb][inb].imag;
            v1 = action.e[ina][inb].imag + action.e[inb][ina].imag;
            v2 = action.e[ina][inb].real - action.e[inb][ina].real;

            // Normalize u
            vsq = v0 * v0 + v1 * v1 + v2 * v2 + v3 * v3;
            z = sqrt((double)vsq);
            v0 = v0/z; v1 = v1/z; v2 = v2/z; v3 = v3/z;
            check = 1.0 - v0 * v0 - v1 * v1 - v2 * v2 - v3 * v3;
            node0_printf("TEST %e ", check);
#endif

            v0 = action.e[ina][ina].real + action.e[inb][inb].real;
            v3 = action.e[ina][ina].imag - action.e[inb][inb].imag;
            v1 = action.e[ina][inb].imag + action.e[inb][ina].imag;
            v2 = action.e[ina][inb].real - action.e[inb][ina].real;
            vsq = v0 * v0 + v1 * v1 + v2 * v2 + v3 * v3;
            z = sqrt((double)vsq);
            norm = 1.0 / sqrt((double)vsq);
            v0 *= norm;
            v1 *= norm;
            v2 *= norm;
            v3 *= norm;
#ifdef DEBUG_PRINT
//            check = 1.0 - v0 * v0 - v1 * v1 - v2 * v2 - v3 * v3;
//            node0_printf("%e\n", check);
//			vsq = v0*v0 + v1*v1 + v2*v2 + v3*v3;
//			z = sqrt((double)vsq );
//			v0 = v0/z; v1 = v1/z; v2 = v2/z; v3 = v3/z;
// test
//node0_printf("v= %e %e %e %e\n",v0,v1,v2,v3);
//node0_printf("z= %e\n",z);
#endif

            /* end norm check--trial SU(2) matrix is a0 + i a(j)sigma(j)*/

            /* test
               node0_printf("v= %.4g %.4g %.4g %.4g\n", v0, v1, v2, v3);
               node0_printf("z= %.4g\n", z);
               */

            // Now begin quasi-heatbath
            // Get four random numbers
            // Add a small increment to avoid log(0)
            xr1 = log((double)(myrand(&(s->site_prn)) + INC));
            xr2 = log((double)(myrand(&(s->site_prn)) + INC));
            xr3 = cos((double)TWOPI * myrand(&(s->site_prn)));
            xr4 = myrand(&(s->site_prn));

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
               rewrite beta/NCOL * re tr(h*v) * z as al*a0
               a0 has prob(a0) = n0 * sqrt(1 - a0**2) * exp(al * a0)
               */
            al = b3 * z;
#ifdef DEBUG_PRINT
            if (lattice[i].x == 1 && lattice[i].y == 2 &&
                lattice[i].z == 0 && lattice[i].t == 1) {
              printf("rand = %.8g %.8g %.8g %.8g and al = %.8g on node %d\n",
                     xr1, xr2, xr3, xr4, al, this_node);
            }
#endif

            /*
               let a0 = 1 - del**2
               get d = del**2
               such that prob2(del) = n1 * del**2 * exp(-al*del**2)
               */

            d = -(xr2  + xr1 * xr3 * xr3) / al;

            /*     monte carlo prob1(del) = n2 * sqrt(1 - 0.5*del**2)
                   then prob(a0) = n3 * prob1(a0)*prob2(a0)
                   */

            // Now beat each site into submission
            nacd = 0;
            if ((1.0 - 0.5 * d) > xr4 * xr4)
              nacd=1;

            // Kennedy--Pendleton algorithm
            if (nacd == 0 && al > 2.0) {
              test=0;
              for (k=0;k<20 && test == 0;k++) {
                kp++;
                /*  get four random numbers
                 *  (add a small increment to prevent taking log(0.)*/
                xr1 = log((double)(myrand(&(s->site_prn)) + INC));
                xr2 = log((double)(myrand(&(s->site_prn)) + INC));
                xr3 = cos((double)TWOPI * myrand(&(s->site_prn)));
                xr4 = myrand(&(s->site_prn));
                d = -(xr2 + xr1 * xr3 * xr3) / al;
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
                xr1 = myrand(&(s->site_prn));
                xr2 = myrand(&(s->site_prn));

                r = xl + xd * xr1;
                a0 = 1.0 + log((double)r) / al;
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
            r2 = fabs(1.0 - a0 * a0);
            r = sqrt(r2);

            /* compute a3 */
            a3 = (2.0*myrand(&(s->site_prn)) - 1.0)*r;

            /* compute a1 and a2 */
            rho = sqrt(fabs(r2 - a3 * a3));

            // xr2 is a random number between 0 and 2pi
            xr2 = TWOPI * myrand(&(s->site_prn));

            a1 = rho * cos((double)xr2);
            a2 = rho * sin((double)xr2);

            // Now update u --> h * u with h = a * v^dag
            h0 = a0 * v0 + a1 * v1 + a2 * v2 + a3 * v3;
            h1 = a1 * v0 - a0 * v1 + a2 * v3 - a3 * v2;
            h2 = a2 * v0 - a0 * v2 + a3 * v1 - a1 * v3;
            h3 = a3 * v0 - a0 * v3 + a1 * v2 - a2 * v1;

            // Elements of SU(2) matrix
            h.e[0][0] = cmplx( h0, h3);
            h.e[0][1] = cmplx( h2, h1);
            h.e[1][0] = cmplx(-h2, h1);
            h.e[1][1] = cmplx( h0,-h3);

            // Update the link
            left_su2_hit_n(&h, ina, inb, &(s->link[dir]));
          }

          // Reunitarize after each SU(2) subgroup sweep
//          reunitarize();

#ifdef DEBUG_PRINT
          // Monitor plaquette after each SU(2) subgroup sweep
          plaquette(&ss_plaq, &st_plaq);
          node0_printf("PLAQ %.8g %.8g %.8g\n",
                       ss_plaq, st_plaq, ss_plaq + st_plaq);
#endif

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
  }
}
// -----------------------------------------------------------------
