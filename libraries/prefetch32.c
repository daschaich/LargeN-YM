// -----------------------------------------------------------------
// Cache prefetching with 8-byte alignment and 32-byte cache line
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

// 8 floats in cache line
// 2 floats per alignment boundary

// matrix: 18 Reals
#define _pftch_M(a) \
  dummy = *((Real *)(a)     ); \
  dummy = *((Real *)(a) +  8); \
  dummy = *((Real *)(a) + 16);

// 4 matrices: 72 Reals
#define _pftch_4M(a) \
  dummy = *((Real *)(a)     ); \
  dummy = *((Real *)(a) +  8); \
  dummy = *((Real *)(a) + 16); \
  dummy = *((Real *)(a) + 24); \
  dummy = *((Real *)(a) + 32); \
  dummy = *((Real *)(a) + 40); \
  dummy = *((Real *)(a) + 48); \
  dummy = *((Real *)(a) + 56); \
  dummy = *((Real *)(a) + 64); \
  dummy = *((Real *)(a) + 70);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void _prefetch_M(matrix *a0) {
  register Real dummy;

  _pftch_M(a0);
}

void _prefetch_4M(matrix *a0) {
  register Real dummy;

  _pftch_4M(a0);
}

#undef _pftch_M
#undef _pftch_4M
// -----------------------------------------------------------------
