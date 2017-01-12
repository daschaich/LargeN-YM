// -----------------------------------------------------------------
// Cache prefetching with 8-byte alignment and 64-byte cache line
// This is the vanilla C version for use when there is no assembly code
// Some compilers will not do these null fetches
// Compiling with the -g option may help with that
// Otherwise you may have to get the compiler-generated assembly code
// and edit by hand

// Fetch addr, addr + cache, ..., addr + size - align
// to get datum of size "size" starting at address "addr"
// with alignment "align" and cache line "cache"

#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/prefetch.h"

// 16 floats in cache line
//  2 floats per alignment boundary

// su3_matrix: 18 Reals
#define _pftch_M(a) \
  dummy = *((Real *)(a)     ); \
  dummy = *((Real *)(a) + 16);

// su3_vector: 6 Reals
#define _pftch_V(a) \
  dummy = *((Real *)(a)     ); \
  dummy = *((Real *)(a) + 4 );

// half_wilson_vector: 12 Reals
#define _pftch_H(a) \
  dummy = *((Real *)(a)     ); \
  dummy = *((Real *)(a) + 10);

// wilson_vector: 24 Reals
#define _pftch_W(a) \
  dummy = *((Real *)(a)     ); \
  dummy = *((Real *)(a) + 16); \
  dummy = *((Real *)(a) + 22);

// 4 su3_vectors: 24 Reals
#define _pftch_4V(a) \
  dummy = *((Real *)(a)     ); \
  dummy = *((Real *)(a) + 16); \
  dummy = *((Real *)(a) + 22);

// 4 su3_matrices: 72 Reals
#define _pftch_4M(a) \
  dummy = *((Real *)(a)     ); \
  dummy = *((Real *)(a) + 16); \
  dummy = *((Real *)(a) + 32); \
  dummy = *((Real *)(a) + 48); \
  dummy = *((Real *)(a) + 64); \
  dummy = *((Real *)(a) + 70);

// 4 wilson_vectors: 96 Reals
#define _pftch_4W(a) \
  dummy = *((Real *)(a)     ); \
  dummy = *((Real *)(a) + 16); \
  dummy = *((Real *)(a) + 32); \
  dummy = *((Real *)(a) + 48); \
  dummy = *((Real *)(a) + 64); \
  dummy = *((Real *)(a) + 80); \
  dummy = *((Real *)(a) + 94);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void _prefetch_M(su3_matrix *a0) {
  register Real dummy;

  _pftch_M(a0);
}

void _prefetch_V(su3_vector *a0) {
  register Real dummy;

  _pftch_V(a0);
}

void _prefetch_W(wilson_vector *a0) {
  register Real dummy;

  _pftch_W(a0);
}

void _prefetch_H(half_wilson_vector *a0) {
  register Real dummy;

  _pftch_H(a0);
}

void _prefetch_VV(su3_vector *a0, su3_vector *a1) {
  register Real dummy;

  _pftch_V(a0);
  _pftch_V(a1);
}

void _prefetch_VVV(su3_vector *a0, su3_vector *a1, su3_vector *a2) {
  register Real dummy;

  _pftch_V(a0);
  _pftch_V(a1);
  _pftch_V(a2);
}

void _prefetch_VVVV(su3_vector *a0, su3_vector *a1, su3_vector *a2,
                    su3_vector *a3) {

  register Real dummy;

  _pftch_V(a0);
  _pftch_V(a1);
  _pftch_V(a2);
  _pftch_V(a3);
}

void _prefetch_VVVVV(su3_vector *a0, su3_vector *a1, su3_vector *a2,
         su3_vector *a3, su3_vector *a4) {
  register Real dummy;

  _pftch_V(a0);
  _pftch_V(a1);
  _pftch_V(a2);
  _pftch_V(a3);
  _pftch_V(a4);
}

void _prefetch_WWW(wilson_vector *a0, wilson_vector *a1, wilson_vector *a2) {
  register Real dummy;

  _pftch_W(a0);
  _pftch_W(a1);
  _pftch_W(a2);
}

void _prefetch_WWWW(wilson_vector *a0, wilson_vector *a1,
        wilson_vector *a2, wilson_vector *a3) {
  register Real dummy;

  _pftch_W(a0);
  _pftch_W(a1);
  _pftch_W(a2);
  _pftch_W(a3);
}

void _prefetch_WWWWW(wilson_vector *a0, wilson_vector *a1, wilson_vector *a2,
                     wilson_vector *a3, wilson_vector *a4) {

  register Real dummy;

  _pftch_W(a0);
  _pftch_W(a1);
  _pftch_W(a2);
  _pftch_W(a3);
  _pftch_W(a4);
}

void _prefetch_4MVVVV(su3_matrix *a0, su3_vector *a1, su3_vector *a2,
                      su3_vector *a3, su3_vector *a4) {

  register Real dummy;

  _pftch_4M(a0);
  _pftch_V(a1);
  _pftch_V(a2);
  _pftch_V(a3);
  _pftch_V(a4);
}

void _prefetch_4MWWWW(su3_matrix *a0, wilson_vector *a1, wilson_vector *a2,
          wilson_vector *a3, wilson_vector *a4) {
  register Real dummy;

  _pftch_4M(a0);
  _pftch_W(a1);
  _pftch_W(a2);
  _pftch_W(a3);
  _pftch_W(a4);
}

void _prefetch_4MV4V(su3_matrix *a0, su3_vector *a1, su3_vector *a2) {
  register Real dummy;

  _pftch_4M(a0);
  _pftch_V(a1);
  _pftch_4V(a2);
}

void _prefetch_4MW4W(su3_matrix *a0, wilson_vector *a1, wilson_vector *a2) {
  register Real dummy;

  _pftch_4M(a0);
  _pftch_W(a1);
  _pftch_4W(a2);
}

#undef _pftch_M
#undef _pftch_V
#undef _pftch_H
#undef _pftch_W
#undef _pftch_4V
#undef _pftch_4M
#undef _pftch_4W
// -----------------------------------------------------------------
