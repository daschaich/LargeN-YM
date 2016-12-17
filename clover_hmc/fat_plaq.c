// -----------------------------------------------------------------
// Do measurements with the smeared NCOLxNCOL link
void fat_plaq() {
  register int i, mu;
  register site *s;
  double ssplaq, stplaq;
  complex plp;

  // Copy the smeared link into s->linkf
  FORALLUPDIR(mu) {
    FORALLSITES(i, s)
      su3mat_copy_f(gauge_field[mu] + i, &(s->linkf[mu]));
  }

  // Measure and print the Polyakov loop and plaquette
  plp = ploop();
  d_plaquette(&ssplaq, &stplaq);
  node0_printf("GFAT %.8g %.8g %.8g %.8g\n",
               (double)plp.real, (double)plp.imag, ssplaq, stplaq);

  // Restore the thin link
  FORALLUPDIR(mu) {
    FORALLSITES(i, s)
      su3mat_copy_f(gauge_field_thin[mu] + i, &(s->linkf[mu]));
  }
}
// -----------------------------------------------------------------
