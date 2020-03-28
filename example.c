#include "geomag.h"
#include <stdio.h>

int main() {
    geomag_init();

    double alt = 0.8; // 0.8 km above WGS84 ellipsoid
    double glat = 40.636; // somewhere around
    double glon = -73.761; // JFK airport
    double time = 2020.0 + (94.0 / 365.0); // day 94 of year 2020
    double dec, dip, ti, gv;
    geomag_calc(alt, glat, glon, time, &dec, &dip, &ti, &gv);
    printf("dec: %.2f\n", dec);
    printf("dip: %.2f\n", dip);
    printf("ti: %.2f\n", ti);
    printf("gv: %.2f\n", gv);

    geomag_interactive();

    geomag_destroy();
}
