/************************ nersc_cksum.c *******************************/
/* MIMD version 7 */

/* Two utilities used in the NERSC archive formate */

/* Compute the low order 32 bits of the unsigned integer sum of the
   float precision real and complex parts of the elements of the gauge
   matrices.
*/

/* Computes the mean global sum of the trace of the gauge links --
   used to aid checking lattice file integrity */

#include "generic_includes.h"

/* commented out as unnecessary [see below], prevent warning [BS 6.07]
static int 
my_big_endian() {
  union  {
    long l;
    char c[sizeof (long)];
  } u;
  u.l = 1;
  return (u.c[sizeof (long) - 1] == 1);
}
*/

u_int32type 
nersc_cksum( void ) {
  u_int32type chksum = 0;

/* all this is commented out so that the checksum will be zero
as initialized above.  [BS 2.07]

  int i,mu,a,b;
  site *s;
  union {
    float       flt;
    u_int32type p32;
  } tmp;
  int big_end = my_big_endian();
  u_int32type p32;

  FORALLSITES(i,s) {
    for(mu=0; mu<4; ++mu) {
      for(a=0; a<NCOL-1; a++) for(b=0; b<NCOL; b++) {
	tmp.flt = s->linkf[mu].e[a][b].real;
	p32 = tmp.p32;
	// if(!big_end)byterevn((int32type *)&p32,1);
	chksum += p32;
	tmp.flt = s->linkf[mu].e[a][b].imag;
	p32 = tmp.p32;
	// if(!big_end)byterevn((int32type *)&p32,1);
	chksum += p32;
      }
    }
  }

  g_uint32sum(&chksum);
*/

  return chksum;

} /* nersc_cksum.c */


void d_linktrsum(double_complex *linktrsum) {
  int i,dir;
  site *s;
  su3_matrix_f *a;

  linktrsum->real = 0.;
  linktrsum->imag = 0.;

  FORALLSITES(i,s){
    FORALLUPDIR(dir){
      a = &s->linkf[dir];
      CSUM(*linktrsum,a->e[0][0]);
      CSUM(*linktrsum,a->e[1][1]);
#if (NCOL>2)
      CSUM(*linktrsum,a->e[2][2]);
#if (NCOL>3)
      CSUM(*linktrsum,a->e[3][3]);
#endif
#endif
    }
  }

  g_dcomplexsum(linktrsum);
  CDIVREAL(*linktrsum,(4*volume),*linktrsum);

} /* d_linktrsum */

