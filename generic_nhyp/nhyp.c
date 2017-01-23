// -----------------------------------------------------------------
// Helper routines for nHYP blocking, specific to SU(2), SU(3), SU(4)
#include "generic_nhyp_includes.h"

#if NCOL == 2
#include "nhyp_SU2.c"
#elif NCOL == 3
#include "nhyp_SU3.c"
#elif NCOL == 4
#include "nhyp_SU4.c"
#else
#include "nhyp_polcof.c"
#endif
// -----------------------------------------------------------------
