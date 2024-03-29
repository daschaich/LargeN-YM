// -----------------------------------------------------------------
// Measure average space--space and space--time plaquettes,
// Uses tempmat for temporary storage
#include "generic_includes.h"

void plaquette(double *ss_plaq, double *st_plaq) {
  register int i, dir, dir2;
  register site *s;
  register matrix *m1, *m4;
  double ss_sum = 0.0, st_sum = 0.0;
  msg_tag *mtag0, *mtag1;
  matrix tmat;

  // We can exploit a symmetry under dir<-->dir2
  for (dir = YUP; dir <= TUP; dir++) {
    for (dir2 = XUP; dir2 < dir; dir2++) {
      // gen_pt[0] is U_b(x+a), gen_pt[1] is U_a(x+b)
      mtag0 = start_gather_site(F_OFFSET(link[dir2]), sizeof(matrix),
                                dir, EVENANDODD, gen_pt[0]);
      mtag1 = start_gather_site(F_OFFSET(link[dir]), sizeof(matrix),
                                dir2, EVENANDODD, gen_pt[1]);

      // tempmat(x) = Udag_b(x) U_a(x)
      FORALLSITES(i, s) {
        m1 = &(s->link[dir]);
        m4 = &(s->link[dir2]);
        mult_an(m4, m1, &(tempmat[i]));
      }
      wait_gather(mtag0);
      wait_gather(mtag1);

      // Compute tr[Udag_a(x+b) Udag_b(x) U_a(x) U_b(x+a)]
      FORALLSITES(i, s) {
        m1 = (matrix *)(gen_pt[0][i]);
        m4 = (matrix *)(gen_pt[1][i]);
        mult_nn(&(tempmat[i]), m1, &tmat);

        if (dir == TUP)
          realtrace_sum(m4, &tmat, &st_sum);
        else
          realtrace_sum(m4, &tmat, &ss_sum);
      }
      cleanup_gather(mtag0);
      cleanup_gather(mtag1);
    }
  }
  g_doublesum(&ss_sum);
  g_doublesum(&st_sum);

  // Average over three plaquettes that involve the temporal link
  // and three that do not
  *ss_plaq = ss_sum / (double)(3.0 * volume);
  *st_plaq = st_sum / (double)(3.0 * volume);
}
// -----------------------------------------------------------------
