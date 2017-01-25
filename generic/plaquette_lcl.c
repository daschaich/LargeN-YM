/* This version mallocs the temporary matrix */

/* Measure the average plaquette of the space-space and
   space-time plaquettes */

/* Added plaq measurement on hypersurfaces; turn on with:          */
/* #define LOCAL_PLAQ */

/* Added measurement of min plaq per config for tuning nhyp as in
Hasenfratz & Knechtli, hep-lat/0103029; turn on with:          */
#define MIN_PLAQ

/* Added local plaq measurement for plotting plaq distribution
 (frep only); turn on with:  */
/* #define ALL_PLAQ */
/*** CAUTION: Do not run ALL_PLAQ with MPI!! ***/

/* The following should be defined self-consistently.
   The value of MY_DIR is printed out on the first time
   the routine is called.                                          */

#ifdef LOCAL_PLAQ
#define MY_X x
#define MY_DIR XUP
#define MY_N nx
static int print_dir=0;
#endif

/* plaq_prll[xx] gives the average of the plaquettes lying in a
   fixed MY_X hypersurface, as a function of the MY_X coordinate.
   Similarly the other plaquettes are given by plaq_prll[xx].      */

#include "generic_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void plaquette_lcl(double *ss_plaq,double *st_plaq) {
  register int i,dir,dir2;
  register site *s;
  register matrix_f *m1, *m4;
  double ss_sum = 0.0, st_sum = 0.0, cur_plaq;
#ifdef MIN_PLAQ
  double min_plaq = NCOL;
#endif
  msg_tag *mtag0,*mtag1;
  matrix_f tmat, *tempmat = malloc(sites_on_node * sizeof(matrix_f));

#ifdef LOCAL_PLAQ
  int xx;
  double *plaq_perp,*plaq_prll;

  plaq_perp = (double *)malloc(MY_N*sizeof(double));
  plaq_prll = (double *)malloc(MY_N*sizeof(double));

  for (xx=0;xx<MY_N;xx++) {
    plaq_perp[xx] = 0.0;
    plaq_prll[xx] = 0.0;
  }
#endif

  if (tempmat == NULL) {
    printf("plaquette: can't malloc tempmat\n");
    fflush(stdout);
    terminate(1);
  }

  // We can exploit a symmetry under dir<-->dir2
  for (dir = YUP; dir <= TUP; dir++) {
    for (dir2 = XUP;dir2 < dir; dir2++) {
      // gen_pt[0] is U_b(x+a), gen_pt[1] is U_a(x+b)
      mtag0 = start_gather_site(F_OFFSET(linkf[dir2]), sizeof(matrix_f),
                                dir, EVENANDODD, gen_pt[0]);
      mtag1 = start_gather_site(F_OFFSET(linkf[dir]), sizeof(matrix_f),
                                dir2, EVENANDODD, gen_pt[1]);

      // tempmat = Udag_b(x) U_a(x)
      FORALLSITES(i,s) {
        m1 = &(s->linkf[dir]);
        m4 = &(s->linkf[dir2]);
        mult_an_f(m4, m1, &tempmat[i]);
      }
      wait_gather(mtag0);
      wait_gather(mtag1);

      // Compute tr[Udag_a(x+b) Udag_b(x) U_a(x) U_b(x+a)]
      FORALLSITES(i,s) {
        m1 = (matrix_f *)(gen_pt[0][i]);
        m4 = (matrix_f *)(gen_pt[1][i]);
        mult_nn_f(&(tempmat[i]), m1, &tmat);
        cur_plaq = (double)realtrace_f(m4, &tmat);
#ifdef MIN_PLAQ
        if (cur_plaq < min_plaq)
          min_plaq = cur_plaq;
#endif
        if (dir == TUP)
          st_sum += cur_plaq;
        else
          ss_sum += cur_plaq;
#ifdef LOCAL_PLAQ
        if (dir == MY_DIR || dir2 == MY_DIR)
          plaq_perp[s->MY_X] += cur_plaq;
        else
          plaq_prll[s->MY_X] += cur_plaq;
#endif
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
  free(tempmat);

#ifdef LOCAL_PLAQ
  for (xx = 0; xx < MY_N; xx++) {
    g_doublesum(&(plaq_perp[xx]));
    g_doublesum(&(plaq_prll[xx]));
  }
  /* normalization */
  for (xx = 0; xx < MY_N; xx++) {
    plaq_perp[xx] *= ((double)(MY_N))/(3.0 * volume);
    plaq_prll[xx] *= ((double)(MY_N))/(3.0 * volume);
  }

  // Print out
  if (this_node==0) {
    if (print_dir==0) {
      printf("LOCAL_PLAQ [0=XUP,..., 3=TUP] dir=%d\n", MY_DIR);
      print_dir=1;
    }
    printf("THIN_PLAQ_PERP");
    for (xx=0;xx<MY_N;xx++)
      printf(" %e", (double)plaq_perp[xx]);
    printf("\n");
    printf("THIN_PLAQ_PRLL");
    for (xx=0;xx<MY_N;xx++)
      printf(" %e", (double)plaq_prll[xx]);
    printf("\n");
  }
#endif

#ifdef MIN_PLAQ
  min_plaq = -min_plaq;
  g_doublemax(&min_plaq);
  min_plaq = -min_plaq;
  node0_printf("MIN_PLAQ_FUND %e\n",min_plaq);
#endif
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Plaquette in fermion irrep
void plaquette_frep_lcl(double *ss_plaq_frep, double *st_plaq_frep) {
  register int i,dir,dir2;
  register site *s;
  register matrix *m1, *m4;
  double ss_sum = 0.0, st_sum = 0.0, cur_plaq;
#ifdef MIN_PLAQ
  double min_plaq = DIMF;
#endif
  msg_tag *mtag0, *mtag1;
  matrix tmat, *tempmat = malloc(sizeof(matrix)*sites_on_node);

#ifdef LOCAL_PLAQ
  int xx;
  double *plaq_perp = malloc(MY_N * sizeof(*plaq_perp));
  double *plaq_prll = malloc(MY_N * sizeof(*plaq_prll));

  for (xx = 0; xx < MY_N; xx++) {
    plaq_perp[xx] = 0.0;
    plaq_prll[xx] = 0.0;
  }
#endif

  if (tempmat == NULL) {
    printf("plaquette: can't malloc tempmat\n");
    fflush(stdout);
    terminate(1);
  }

  for (dir = YUP; dir <= TUP; dir++) {
    for (dir2 = XUP; dir2 < dir; dir2++) {

      mtag0 = start_gather_site(F_OFFSET(link[dir2]), sizeof(matrix),
          dir, EVENANDODD, gen_pt[0]);
      mtag1 = start_gather_site(F_OFFSET(link[dir]), sizeof(matrix),
          dir2, EVENANDODD, gen_pt[1]);

      FORALLSITES(i,s) {
        m1 = &(s->link[dir]);
        m4 = &(s->link[dir2]);
        mult_an(m4,m1,&tempmat[i]);
      }
      wait_gather(mtag0);
      wait_gather(mtag1);

      if (dir==TUP) {
        FORALLSITES(i, s) {
          mult_nn(&(tempmat[i]), (matrix *)(gen_pt[0][i]),
              &tmat);
          cur_plaq = (double)
            realtrace((matrix *)(gen_pt[1][i]),&tmat);
#ifdef MIN_PLAQ
          if (cur_plaq<min_plaq) min_plaq=cur_plaq;
#endif
#ifdef ALL_PLAQ
          printf("ALL_PLAQ %d %d %d %d %d %d %e\n",
              s->x,s->y,s->z,s->t,dir,dir2,cur_plaq);
#endif
          st_sum += cur_plaq;
#ifdef LOCAL_PLAQ
          if (dir==MY_DIR || dir2==MY_DIR) {
            plaq_perp[s->MY_X] += cur_plaq;
          }
          else{
            plaq_prll[s->MY_X] += cur_plaq;
          }
#endif
        }
      }
      else {
        FORALLSITES(i,s) {
          mult_nn(&(tempmat[i]), (matrix *)(gen_pt[0][i]),
              &tmat);
          cur_plaq = (double)
            realtrace((matrix *)(gen_pt[1][i]),&tmat);
#ifdef MIN_PLAQ
          if (cur_plaq<min_plaq) min_plaq=cur_plaq;
#endif
#ifdef ALL_PLAQ
          printf("ALL_PLAQ %d %d %d %d %d %d %e\n",
              s->x,s->y,s->z,s->t,dir,dir2,cur_plaq);
#endif
          ss_sum += cur_plaq;
#ifdef LOCAL_PLAQ
          if (dir==MY_DIR || dir2==MY_DIR)
            plaq_perp[s->MY_X] += cur_plaq;
          else
            plaq_prll[s->MY_X] += cur_plaq;
#endif
        }
      }
      cleanup_gather(mtag0);
      cleanup_gather(mtag1);
    }
  }
  g_doublesum(&ss_sum);
  g_doublesum(&st_sum);
  // Average over three plaquettes that involve the temporal link
  // and three that do not
  *ss_plaq_frep = ss_sum / (double)(3.0 * volume);
  *st_plaq_frep = st_sum / (double)(3.0 * volume);

  free(tempmat);

#ifdef LOCAL_PLAQ
  for (xx = 0; xx < MY_N; xx++) {
    g_doublesum(&(plaq_perp[xx]));
    g_doublesum(&(plaq_prll[xx]));
  }

  // normalization
  for (xx = 0; xx < MY_N; xx++) {
    plaq_perp[xx] *= (double)(MY_N) / (3.0 * volume);
    plaq_prll[xx] *= (double)(MY_N) / (3.0 * volume);
  }

  // Print out
  if (this_node==0) {
    /*    printf("LOCAL_PLAQ dir=%d\n", MY_DIR); */
    printf("FAT_PLAQ_PERP");
    for (xx = 0; xx < MY_N; xx++)
      printf(" %e", (double)plaq_perp[xx]);
    printf("\n");
    printf("FAT_PLAQ_PRLL");
    for (xx = 0; xx < MY_N; xx++)
      printf(" %e", (double)plaq_prll[xx]);
    printf("\n");
  }
#endif /* LOCAL_PLAQ */

#ifdef MIN_PLAQ
  min_plaq = -min_plaq;
  g_doublemax(&min_plaq);
  min_plaq = -min_plaq;
  node0_printf("MIN_PLAQ_FERM %e\n",min_plaq);
#endif
}
// -----------------------------------------------------------------
