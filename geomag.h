
void geomag_init();
void geomag_destroy();

void geomag_calc(double alt, double glat, double glon, double time, double *dec, double *dip, double *ti, double *gv);
void geomag_interactive();
