/**************  r_su2_hit_a.c  (in su3.a) **********************
*									*
*  right multiply an matrix by the adjoint of an su2 matrix 	*
*/

#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"

void right_su2_hit_a(su2_matrix *u,int p,int q,matrix *link)
{
  /* link <-  link * u adj */
  /* The 0 column of u-adjoint matches column p of the SU(3) matrix */
  /* The 1 column of u-adjoint matches column q of the SU(3) matrix */
  /* C. DeTar 18 Oct 1990 */

  register int m;

  for (m = 0; m < DIMF; m++)
    mult_su2_mat_vec_elem_a(u, &(link->e[m][p]), &(link->e[m][q]));

} /* r_su2_hit_a.c */

