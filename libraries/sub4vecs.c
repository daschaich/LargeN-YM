// -----------------------------------------------------------------
// Subtract four su3 vectors
// a <-- a - b1 - b2 - b3 - b4
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void sub_four_vecs(vector *a, vector *b1, vector *b2,
                   vector *b3, vector *b4) {
#ifndef FAST
  register int i;
  for (i = 0; i < DIMF; i++) {
    CDIF(a->c[i], b1->c[i]);
    CDIF(a->c[i], b2->c[i]);
    CDIF(a->c[i], b3->c[i]);
    CDIF(a->c[i], b4->c[i]);
  }
#else  // FAST version for NCOL = DIMF = 3
  CDIF(a->c[0], b1->c[0]);
  CDIF(a->c[1], b1->c[1]);
  CDIF(a->c[2], b1->c[2]);
  CDIF(a->c[0], b2->c[0]);
  CDIF(a->c[1], b2->c[1]);
  CDIF(a->c[2], b2->c[2]);
  CDIF(a->c[0], b3->c[0]);
  CDIF(a->c[1], b3->c[1]);
  CDIF(a->c[2], b3->c[2]);
  CDIF(a->c[0], b4->c[0]);
  CDIF(a->c[1], b4->c[1]);
  CDIF(a->c[2], b4->c[2]);
#endif
}
// -----------------------------------------------------------------
