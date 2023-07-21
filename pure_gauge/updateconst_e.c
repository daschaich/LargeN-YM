// -----------------------------------------------------------------
// Update lattice
// Choose between hybrid Monte Carlo trajectories
// or over-relaxed quasi-heatbath (qhb) sweeps
#include "pg_includes.h"

void updateconst_e(double E_min) {
#ifdef LLR
  // Check unitarity before doing anything
  check_unitarity();

#ifdef HMC
  // Do HMC updates if specified
  //update_hmc(E_min);
  //if(this_node==0){
    //printf("Node nr %d is doing heatbath\n", this_node);
  //monteconst_e();
  int i;
  for(i=0;i<numnodes();i++){
    frankensteinconst_e(E_min,i);
    //monte();
  }
  
  
  //printf("Node nr %d is doing heatbath\n", this_node);
  //}
  //else{
  //  printf("Node nr %d is doing overrelax\n", this_node);
  //  relax();
  //}
#else
  // Otherwise do over-relaxation and quasi-heatbath steps
  relax();
  monteconst_e();
#endif

  // Reunitarize the gauge field
  reunitarize();
#endif
}
// -----------------------------------------------------------------
