// -----------------------------------------------------------------
// Kennedy--Pendleton quasi-heatbath (qhb) on SU(2) subgroups and Microcanonical over-relaxation by doing successive SU(2) gauge hits
#include "pg_includes.h"
#define INC 1.0e-10
//#define DEBUG_PRINT

void frankensteinconst_e(double E_min,int node_counter) {

//heatbath & overrelax
  register int dir, i;
  register site *s;
  int istep, Nhit, subgrp, ina, inb, parity, count;
  int k, kp, cr, nacd, test, index_a[N_OFFDIAG], index_b[N_OFFDIAG];
  int this_accept = 0, this_reject = 0;
  Real xr1, xr2, xr3, xr4;
  Real a0 = 0, a1, a2, a3;
  Real v0, v1, v2, v3, vsq;
  Real h0, h1, h2, h3;
  Real r, r2, rho, z, norm;
  Real al, d, xl, xd, b3 = beta * a * one_ov_N;
  double E, ss_plaq, st_plaq;
  double Etest, Etestb4, Etestafter,Etrb4,Etrafter, Edif;
  su2_matrix h;
  matrix actmat;
  
//remaining overrelax var
  Real asq;
  su2_matrix u;
  matrix action; 
  
  
  
  
#ifdef DEBUG_PRINT
  double rate, check;
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
  
  E = gauge_action();
  Etest = E;

  // Loop over sweeps
  for (istep = 0; istep < qhb_steps; istep++) {
    for (parity = ODD; parity <= EVEN; parity++) {
      FORALLUPDIR(dir) {
        // Compute the gauge force (updating s->staple)
        dsdu_qhb(dir, parity);

        // Now for the qhb&or updating, looping over SU(2) subgroups
        for (subgrp = 0; subgrp < Nhit; subgrp++) {
          kp = 0;
          cr = 0;

          // Pick out this SU(2) subgroup
          ina = index_a[subgrp];
          inb = index_b[subgrp];
         
         //node0_printf("Energy b4 OR and HB %.8g\n", E);
//overrelax
        if(this_node!=node_counter){
          //printf("Node nr %d is doing overrelax \n", this_node);
          FORSOMEPARITY(i, s, parity) {
            mult_na(&(s->link[dir]), &(s->staple), &action);

            // Decompose the action into SU(2) subgroups
            // using Pauli matrix expansion
            // The SU(2) hit matrix is represented as
            //    a0 + i * Sum j (sigma j * aj)
            mult_na(&(s->link[dir]), &(s->staple), &action);
#ifdef DEBUG_PRINT
            a0 = action.e[ina][ina].real + action.e[inb][inb].real;
            a3 = action.e[ina][ina].imag - action.e[inb][inb].imag;
            a1 = action.e[ina][inb].imag + action.e[inb][ina].imag;
            a2 = action.e[ina][inb].real - action.e[inb][ina].real;

            // Normalize and complex conjugate u
            asq = a0 * a0 + a1 * a1 + a2 * a2 + a3 * a3;
            r = sqrt((double)asq);
            a0 = a0/r; a1 = -a1/r; a2 = -a2/r; a3 = -a3/r;
            test = 1.0 - a0 * a0 - a1 * a1 - a2 * a2 - a3 * a3;
            node0_printf("TEST %e ", test);
#endif

            a0 = action.e[ina][ina].real + action.e[inb][inb].real;
            a3 = action.e[ina][ina].imag - action.e[inb][inb].imag;
            a1 = action.e[ina][inb].imag + action.e[inb][ina].imag;
            a2 = action.e[ina][inb].real - action.e[inb][ina].real;
            asq = a0 * a0 + a1 * a1 + a2 * a2 + a3 * a3;
            norm = 1.0 / sqrt((double)asq);
            a0 *=  norm;
            a1 *= -norm;
            a2 *= -norm;
            a3 *= -norm;

            // Elements of SU(2) matrix
            u.e[0][0] = cmplx( a0, a3);
            u.e[0][1] = cmplx( a2, a1);
            u.e[1][0] = cmplx(-a2, a1);
            u.e[1][1] = cmplx( a0,-a3);

            // Do SU(2) hit on all links twice (to overrelax)
            left_su2_hit_n(&u, ina, inb, &(s->link[dir]));
            left_su2_hit_n(&u, ina, inb, &(s->link[dir]));
          }
          }
          //E = gauge_action();
          //node0_printf("Energy after OR %.8g\n", E);
//heatbath
          if(this_node==node_counter){
          //printf("Node nr %d is doing heatbath\n", this_node);
          FORSOMEPARITY(i, s, parity) {
            // Save starting links
            mat_copy(&(s->link[dir]), &(s->tempmat));

            // Decompose the action into SU(2) subgroups
            // using Pauli matrix expansion
            // The SU(2) hit matrix is represented as
            //   v0 + i * Sum j (sigma j * vj)
            mult_na(&(s->link[dir]), &(s->staple), &actmat);
#ifdef DEBUG_PRINT
            v0 = actmat.e[ina][ina].real + actmat.e[inb][inb].real;
            v3 = actmat.e[ina][ina].imag - actmat.e[inb][inb].imag;
            v1 = actmat.e[ina][inb].imag + actmat.e[inb][ina].imag;
            v2 = actmat.e[ina][inb].real - actmat.e[inb][ina].real;

            // Normalize u
            vsq = v0 * v0 + v1 * v1 + v2 * v2 + v3 * v3;
            z = sqrt((double)vsq);
            v0 = v0/z; v1 = v1/z; v2 = v2/z; v3 = v3/z;
            check = 1.0 - v0 * v0 - v1 * v1 - v2 * v2 - v3 * v3;
            node0_printf("TEST %e ", check);
#endif

            v0 = actmat.e[ina][ina].real + actmat.e[inb][inb].real;
            v3 = actmat.e[ina][ina].imag - actmat.e[inb][inb].imag;
            v1 = actmat.e[ina][inb].imag + actmat.e[inb][ina].imag;
            v2 = actmat.e[ina][inb].real - actmat.e[inb][ina].real;
            vsq = v0 * v0 + v1 * v1 + v2 * v2 + v3 * v3;
            z = sqrt((double)vsq);
            norm = 1.0 / sqrt((double)vsq);
            v0 *= norm;
            v1 *= norm;
            v2 *= norm;
            v3 *= norm;

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
               rewrite beta/3 * re tr(h*v) * z as al*a0
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
            Etrb4=-realtrace(&(s->link[dir]),&(s->staple));
            left_su2_hit_n(&h, ina, inb, &(s->link[dir]));
            Etrafter=-realtrace(&(s->link[dir]),&(s->staple));
            Edif=Etrafter-Etrb4;
            Etest=Etest+one_ov_N*Edif;
#ifdef DEBUG_PRINT
          printf("PLAQ %.8g %.8g %.8g\n",
                       ss_plaq, st_plaq, ss_plaq + st_plaq);
#endif

          if (Etest < E_min || Etest > E_min+delta) {
            this_reject++;
            //printf("Energy %.8g leaves [%.8g, %.8g]\n", Etest, E_min, E_min+delta);
            Etest=Etest-one_ov_N*Edif;
#ifdef DEBUG_PRINT
            printf("Reject subgrp %d for parity %d: ", subgrp, parity);
            printf("Energy %.8g leaves [%.8g, %.8g]\n", Etest, E_min, E_min+delta);
#endif
            
            mat_copy(&(s->tempmat), &(s->link[dir]));
          }
          else {
            this_accept++;
            //node0_printf("Accept new energy %.8g\n", E);
            //printf("Accept new energy %.8g\n", Etest);
#ifdef DEBUG_PRINT
            node0_printf("Accept new energy %.8g\n", E);
            printf("Accept new energy %.8g\n", Etest);
#endif
          }
            
          }
          }
          // Reunitarize after each SU(2) subgroup sweep
//          reunitarize();

          // If we have exited the energy interval, restore starting links
          //plaquette(&ss_plaq, &st_plaq);
          
          //E = gauge_action();
          //node0_printf("Energy after HB %.8g\n", E);
          //node0_printf("Etest after HB %.8g\n", Etest);
          
          // Monitor plaquette after each SU(2) subgroup sweep
/*          if(this_node==node_counter){
#ifdef DEBUG_PRINT
          printf("PLAQ %.8g %.8g %.8g\n",
                       ss_plaq, st_plaq, ss_plaq + st_plaq);
#endif
          if (E < E_min || E > E_min+delta) {
            this_reject++;
            printf("Energy %.8g leaves [%.8g, %.8g]\n", E, E_min, E_min+delta);
            Etest=E;
#ifdef DEBUG_PRINT
            printf("Reject subgrp %d for parity %d: ", subgrp, parity);
            printf("Energy %.8g leaves [%.8g, %.8g]\n", E, E_min, E_min+delta);
#endif
            FORSOMEPARITY(i, s, parity)
              mat_copy(&(s->tempmat), &(s->link[dir]));
          }
          else {
            this_accept++;
            //node0_printf("Accept new energy %.8g\n", E);
            printf("Accept new energy %.8g\n", E);
#ifdef DEBUG_PRINT
            //node0_printf("Accept new energy %.8g\n", E);
            printf("Accept new energy %.8g\n", E);
#endif
          }
          }
*/          
        }
      }
    }
  }
  // Update overall acceptance
  accept += this_accept;
  reject += this_reject;
  node0_printf("Energy after Updates %.8g\n", Etest);

#ifdef DEBUG_PRINT
  rate = (double)this_accept / ((double)(this_accept + this_reject));
  node0_printf("Acceptance %d of %d = %.4g\n",
               this_accept, this_accept + this_reject, rate);
  
#endif
}
// -----------------------------------------------------------------
