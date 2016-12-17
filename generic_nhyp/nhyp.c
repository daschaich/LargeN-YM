/* Reference:
* Hypercubic smeared links for dynamical fermions.
* By Anna Hasenfratz, Roland Hoffmann, Stefan Schaefer.
* JHEP 0705:029,2007. [hep-lat/0702028]
*/

/* this file contains helper routines for the stout blocking
*  which compute certain functions of a hermitian 3x3 Matrix Q
*  They follow closely the description in Morningstar & Peardon
*/

#include "generic_nhyp_includes.h"

/* These helper routines are specific to SU(2,3,4) */

void dump_math_f(su3_matrix_f *Q);

#ifdef NHYP_DEBUG
#define DUMP_STUFF               \
    printf("\nOmega is:\n\n");   \
    dump_math_f(Omega);          \
    printf("\nQ is:\n\n");       \
    dump_math_f(Q);              \
    printf("\n");                \
    fflush(stdout);
#endif

#if NCOL == 2
#include "nhyp_SU2.c"
#elif NCOL == 3
#include "nhyp_SU3.c"
#elif NCOL == 4
#include "nhyp_SU4.c"
#else
#include "nhyp_polcof.c"
#endif




/* dumping an NCOLxNCOL matrix in a format closest to Mathematica.
To complete the conversion, run sed with the instruction

s/e\([+-][0-9]*\)[)]/)*10^(\1)/g

*/

#ifdef NHYP_DEBUG
void dump_math_f(su3_matrix_f *Q){
    register int i,j;
    printf("{");
    for(i=0;i<NCOL;i++){
        printf("{");
        for(j=0;j<NCOL;j++){
            printf("(%e)+I*(%e)", Q->e[i][j].real,Q->e[i][j].imag);
            if(j<NCOL-1){
                printf(",\t");
            }
            else{
                if(i<NCOL-1){
                     printf("},\n");
                }
                else{
                     printf("}}\n");
                }
            }
        }
    }
}
#endif
