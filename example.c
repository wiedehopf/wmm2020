#include "geomag.h"
#include <stdio.h>

int main() {
    geomag_init();

    double alt = 0;
    double glat = 40;
    double glon = -77;
    double time = 2020;
    double dec, dip, ti, gv;
    geomag_calc(alt, glat, glon, time, &dec, &dip, &ti, &gv);
    printf("dec: %.2f\n", dec);
    printf("dip: %.2f\n", dip);
    printf("ti: %.2f\n", ti);
    printf("gv: %.2f\n", gv);

    geomag_interactive();

    geomag_destroy();
}
