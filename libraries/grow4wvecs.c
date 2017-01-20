// a <-- b1 + b2 + b3 + b4
// a <-- a + b1 + b2 + b3 + b4

// B1 is expanded using gamma_x, B2 using gamma_y, etc.
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/dirs.h"

void grow_add_four_wvecs(wilson_vector *a, half_wilson_vector *b1,
                         half_wilson_vector *b2, half_wilson_vector *b3,
                         half_wilson_vector *b4, int sign) {
#ifndef FAST
    wp_grow(b1, a, XUP, sign);
    wp_grow_sum(b2, a, YUP, sign);
    wp_grow_sum(b3, a, ZUP, sign);
    wp_grow_sum(b4, a, TUP, sign);

#else  // Inline wp_grow and wp_grow_sum
  int i;

  /* wp_grow(b1, a, XUP, sign); */
  /* case XUP: */
  if (sign == PLUS) {
    for (i = 0; i < DIMF; i++) {
      a->d[0].c[i] = b1->h[0].c[i];
      a->d[1].c[i] = b1->h[1].c[i];
      CMUL_MINUS_I(b1->h[0].c[i], a->d[3].c[i]);
      CMUL_MINUS_I(b1->h[1].c[i], a->d[2].c[i]);
    }
  }
  else {
    /* case XDOWN: */
    for (i = 0; i < DIMF; i++) {
      a->d[0].c[i] = b1->h[0].c[i];
      a->d[1].c[i] = b1->h[1].c[i];
      CMUL_I(b1->h[0].c[i], a->d[3].c[i]);
      CMUL_I(b1->h[1].c[i], a->d[2].c[i]);
    }
  }

  /* wp_grow_sum(b2, a, YUP, sign); */
  if (sign == PLUS) {
    /* case YUP: */
    for (i = 0; i < DIMF; i++) {
      CSUM(a->d[0].c[i], b2->h[0].c[i]);
      CSUM(a->d[1].c[i], b2->h[1].c[i]);
      CSUM(a->d[2].c[i], b2->h[1].c[i]);
      CDIF(a->d[3].c[i], b2->h[0].c[i]);
    }
  }
  else {
    /* case YDOWN: */
    for (i = 0; i < DIMF; i++) {
      CSUM(a->d[0].c[i], b2->h[0].c[i]);
      CSUM(a->d[1].c[i], b2->h[1].c[i]);
      CDIF(a->d[2].c[i], b2->h[1].c[i]);
      CSUM(a->d[3].c[i], b2->h[0].c[i]);
    }
  }

  /* wp_grow_sum(b3, a, ZUP, sign); */
  if (sign == PLUS) {
    /* case ZUP: */
    for (i = 0; i < DIMF; i++) {
      CSUM(a->d[0].c[i], b3->h[0].c[i]);
      CSUM(a->d[1].c[i], b3->h[1].c[i]);
      CSUM_TMI(a->d[2].c[i], b3->h[0].c[i]);
      CSUM_TPI(a->d[3].c[i], b3->h[1].c[i]);
    }
  }
  else {
    /* case ZDOWN:*/
    for (i = 0; i < DIMF; i++) {
      CSUM(a->d[0].c[i], b3->h[0].c[i]);
      CSUM(a->d[1].c[i], b3->h[1].c[i]);
      CSUM_TPI(a->d[2].c[i], b3->h[0].c[i]);
      CSUM_TMI(a->d[3].c[i], b3->h[1].c[i]);
    }
  }

  /* wp_grow_sum(b4, a, TUP, sign); */
  if (sign == PLUS) {
    /* case TUP: */
    for (i = 0; i < DIMF; i++) {
      CSUM(a->d[0].c[i], b4->h[0].c[i]);
      CSUM(a->d[1].c[i], b4->h[1].c[i]);
      CSUM(a->d[2].c[i], b4->h[0].c[i]);
      CSUM(a->d[3].c[i], b4->h[1].c[i]);
    }
  }
  else {
    /* case TDOWN: */
    for (i = 0; i < DIMF; i++) {
      CSUM(a->d[0].c[i], b4->h[0].c[i]);
      CSUM(a->d[1].c[i], b4->h[1].c[i]);
      CDIF(a->d[2].c[i], b4->h[0].c[i]);
      CDIF(a->d[3].c[i], b4->h[1].c[i]);
    }
  }
#endif
}
void grow_sum_four_wvecs(wilson_vector *a, half_wilson_vector *b1,
                         half_wilson_vector *b2, half_wilson_vector *b3,
                         half_wilson_vector *b4, int sign) {
#ifndef FAST
    wp_grow_sum(b1, a, XUP, sign);
    wp_grow_sum(b2, a, YUP, sign);
    wp_grow_sum(b3, a, ZUP, sign);
    wp_grow_sum(b4, a, TUP, sign);

#else  // Inline wp_grow and wp_grow_sum
  int i;

  /* wp_grow_add(b1, a, XUP, sign); */
  /* case XUP: */
  if (sign == PLUS) {
    for (i = 0; i < DIMF; i++) {
      CSUM(a->d[0].c[i], b1->h[0].c[i]);
      CSUM(a->d[1].c[i], b1->h[1].c[i]);
      CSUM_TMI(a->d[2].c[i], b1->h[1].c[i]);
      CSUM_TMI(a->d[3].c[i], b1->h[0].c[i]);
    }
  }
  else {
    /* case XDOWN: */
    for (i = 0; i < DIMF; i++) {
      CSUM(a->d[0].c[i], b1->h[0].c[i]);
      CSUM(a->d[1].c[i], b1->h[1].c[i]);
      CSUM_TPI(a->d[2].c[i], b1->h[1].c[i]);
      CSUM_TPI(a->d[3].c[i], b1->h[0].c[i]);
    }
  }

  /* wp_grow_sum(b2, a, YUP, sign); */
  if (sign == PLUS) {
    /* case YUP: */
    for (i = 0; i < DIMF; i++) {
      CSUM(a->d[0].c[i], b2->h[0].c[i]);
      CSUM(a->d[1].c[i], b2->h[1].c[i]);
      CSUM(a->d[2].c[i], b2->h[1].c[i]);
      CDIF(a->d[3].c[i], b2->h[0].c[i]);
    }
  }
  else {
    /* case YDOWN: */
    for (i = 0; i < DIMF; i++) {
      CSUM(a->d[0].c[i], b2->h[0].c[i]);
      CSUM(a->d[1].c[i], b2->h[1].c[i]);
      CDIF(a->d[2].c[i], b2->h[1].c[i]);
      CSUM(a->d[3].c[i], b2->h[0].c[i]);
    }
  }

  /* wp_grow_sum(b3, a, ZUP, sign); */
  if (sign == PLUS) {
    /* case ZUP: */
    for (i = 0; i < DIMF; i++) {
      CSUM(a->d[0].c[i], b3->h[0].c[i]);
      CSUM(a->d[1].c[i], b3->h[1].c[i]);
      CSUM_TMI(a->d[2].c[i], b3->h[0].c[i]);
      CSUM_TPI(a->d[3].c[i], b3->h[1].c[i]);
    }
  }
  else {
    /* case ZDOWN:*/
    for (i = 0; i < DIMF; i++) {
      CSUM(a->d[0].c[i], b3->h[0].c[i]);
      CSUM(a->d[1].c[i], b3->h[1].c[i]);
      CSUM_TPI(a->d[2].c[i], b3->h[0].c[i]);
      CSUM_TMI(a->d[3].c[i], b3->h[1].c[i]);
    }
  }

  /* wp_grow_sum(b4, a, TUP, sign); */
  if (sign == PLUS) {
    /* case TUP: */
    for (i = 0; i < DIMF; i++) {
      CSUM(a->d[0].c[i], b4->h[0].c[i]);
      CSUM(a->d[1].c[i], b4->h[1].c[i]);
      CSUM(a->d[2].c[i], b4->h[0].c[i]);
      CSUM(a->d[3].c[i], b4->h[1].c[i]);
    }
  }
  else {
    /* case TDOWN: */
    for (i = 0; i < DIMF; i++) {
      CSUM(a->d[0].c[i], b4->h[0].c[i]);
      CSUM(a->d[1].c[i], b4->h[1].c[i]);
      CDIF(a->d[2].c[i], b4->h[0].c[i]);
      CDIF(a->d[3].c[i], b4->h[1].c[i]);
    }
  }
#endif
}
