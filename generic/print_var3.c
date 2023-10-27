// -----------------------------------------------------------------
// Print all complex traces of tempmat on timeslice t=0
// Current specialized to Polyakov loop, but easy to generalize
#include "generic_includes.h"

void print_var3() {
#ifdef LOCALPOLY
  int currentnode = 0, newnode;
  int i, x, y, z, t = 0, node0 = 0;
  complex toprint;

  g_sync();
  for (z = 0; z < nx; z++) {
    for (y = 0; y < ny; y++) {
      for (x = 0; x < nx; x++) {
        // The data to print
        i = node_index(x, y, z, t);
        toprint = trace(&(tempmat[i]));

        // Check whether currentnode has changed as we loop over (x, y, z)
        newnode = node_number(x, y, z, t);
        if (newnode != currentnode) {
          g_sync();
          currentnode = newnode;
        }

        // this_node=mynode() should be set up by each application's setup.c
        if (this_node == 0) {
          if (currentnode != 0)             // Gather from currentnode
            get_field((char *)&toprint, sizeof(complex), currentnode);

          // Now we can print
          if ((printf("POLYDIST %d %d %d %.8g %.8g\n", x, y, z,
                      toprint.real, toprint.imag) == EOF)) {
            printf("print_var3: Write error\n");
            terminate(1);
          }
        }
        else {  // this_node is not node0
          if (this_node == currentnode)     // Gather to node0
            send_field((char *)&toprint, sizeof(complex), node0);
        }
      }
    }
  }
  g_sync();
  if (this_node == 0)
    fflush(stdout);
#endif
}
// -----------------------------------------------------------------
