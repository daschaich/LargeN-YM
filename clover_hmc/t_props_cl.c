/* Measure propagators for measuring m_quark by method of
   Iwasaki et al.

   Only works with LU-preconditioned matrix
*/
#include "cl_dyn_includes.h"

int t_props_cl( ) {
  register site *st;
  register int i, spin, iters = 0;
  Real *piprop = malloc(nt * sizeof(*piprop));
  Real *aprop = malloc(nt * sizeof(*aprop));
  complex cc;
  wilson_vector tvec;

  /* spin zero color zero for now */
  spin=0;

  /* Make source directly in chi[0] */
  FORALLSITES(i, st)
    clear_wvec(&(st->chi[0]));

  /* Wall source */
  /* Even sites only for now */
  FOREVENSITES(i, st){
    if (st->t == 0)
      st->chi[0].d[spin].c[0].real = 1.0;
  }

  // Load inversion control structure
  qic.start_flag = 0;   // Use zero initial guess for psi[0]

  // Load Dirac matrix parameters, including temporaries
  qic.wv1 = F_OFFSET(tmp);
  qic.wv2 = F_OFFSET(mp);

  // Invert M, result in psi[0]
  iters += wilson_invert(F_OFFSET(chi[0]), F_OFFSET(psi[0]), F_OFFSET(r),
      cgilu_cl, &qic, (void *)&dcp);

  /* Sum pion propagator and gamma_3 propagator over z slices  */
  for(i=0;i<nt;i++) piprop[i] = aprop[i] = 0.0;
  FORALLSITES(i,st){
    piprop[st->t] += magsq_wvec( &(st->psi[0]) );
    mult_by_gamma( &(st->psi[0]), &tvec, TUP );
    cc = wvec_dot( &(st->psi[0]), &tvec );
    aprop[st->t] += cc.real;
  }
  for(i=0;i<nt;i++){
    g_floatsum( &(piprop[i]) );
    g_floatsum( &(aprop[i]) );
  }

  /* Dump results */
  if(this_node==0)for(i=0;i<nt;i++)
    printf("TPROPS %d %d %.7e %.7e\n",spin,i,(double)piprop[i],
        (double)aprop[i]);

  free(piprop); free(aprop);
  free_clov();
  return(iters);
}
