/*  uses staple.e[0][0] to accumulate, and ordinary gathers*/

#include "cl_dyn_includes.h"

void gauge_action(double *result) {
  /* double sum_tr; */
  register int i, k;
  register site *s;
  int rep, dirs[10], sign[10], length, ln, iloop;
  double g_action = 0.0, action, act2, total_action;
  complex trace;

  for (iloop = 0; iloop < nloop; iloop++) {
    length = loop_length[iloop];
    /* loop over rotations and reflections */
    for (ln = 0; ln < loop_num[iloop]; ln++) {
      /* set up dirs and sign */
      for (k = 0;k<length;k++) {
        if (loop_table[iloop][ln][k] < 4) {
          sign[k] = 1;
          dirs[k] = (loop_table[iloop][ln][k]) % 4;
        }
        else {
          sign[k] = -1;
          dirs[k] = (7 - loop_table[iloop][ln][k]) % 4;
        }
      }

      path(dirs, sign, length);   // Puts result in tempmatf

      FORALLSITES(i, s) {
        trace = trace_f(&(tempmatf[i]));
        action = (Real)NCOL - (double)trace.real;
        total_action = (double)loop_coeff[iloop][0]*action;
        act2 = action;
        for (rep = 1; rep < nreps; rep++) {
          act2 *= action;
          total_action += (double)loop_coeff[iloop][rep]*act2;
        }
        g_action  += total_action;
      }
    }
  }
  g_doublesum(&g_action);

  *result = g_action;
}
