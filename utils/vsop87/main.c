/*
    This file is part of libamspace.

    libamspace is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    libamspace is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <libamspace.h>
#include <vsop87.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

static const char *s_datadir;

static void printUsage(void)
{
    fprintf(stderr, "Usage: vsop87 <datadir>\n");
}
static inline void printVect(DT v[3])
{
    printf("%0.9f %0.9f %0.9f\n", v[0], v[1], v[2]);
}
static void printBody(vsop87_body_et body, char *name, double T, DT v_e[3],
                                                                double e)
{
    DT v[3], v_g[3], v_ec[3], v_eq[3];
    sexangle_st ang;
    if (vsop87_getCoords(v, s_datadir, VSOP87_HELIO_RECT_DATE, body, T)!=0)
        return;
    printf("%s\theliocentric:\t", name); printVect(v);
    printf("\tgeocentric:\t"); VECT3_DIFF(v_g, v, v_e); printVect(v_g);

    printf("\tecliptic geo:\t");
    amsp_coord_rect2ecliptic(v_ec, v_g);
    amsp_angle_deg2sex(&ang, v_ec[0]*CRG); amsp_angle_sexprint(&ang);
    printf(" ");
    amsp_angle_deg2sex(&ang, v_ec[1]*CRG); amsp_angle_sexprint(&ang);
    printf(" ");
    printf("%0.9lf\n", v_ec[2]);

    printf("\tequatorial geo:\t");
    amsp_coord_ecliptic2equatorial(v_eq, v_ec, e);
    amsp_angle_deg2sex(&ang, v_eq[0]*CRG); amsp_angle_sexprint(&ang);
    printf(" ");
    amsp_angle_deg2sex(&ang, v_eq[1]*CRG); amsp_angle_sexprint(&ang);
    printf(" ");
    printf("%0.9lf\n", v_eq[2]);
}

int main(const int argc, char * const *argv)
{
    dtm_st dtm;
    time_t t;
    DT T, e;
    DT v_e[3], v_g[3];
    sexangle_st ang;

    if (argc!=2) {
        printUsage();
        return EXIT_FAILURE;
    }
    s_datadir = argv[1];

    t = time(NULL);
    gmtime_r(&t, &dtm.tm);
    dtm.sec = 0.;
    T = amsp_date_tm2julian(&dtm);

    e = amsp_coord_earthEclipticAtDate(T);
    printf("Earth ecliptic: %f ", e);
    amsp_angle_deg2sex(&ang, e*CRG);
    amsp_angle_sexprint(&ang); printf("\n");

    if (vsop87_getCoords(v_e, s_datadir, VSOP87_HELIO_RECT_DATE, BODY_EARTH, T)!=0)
        return EXIT_FAILURE;
    printf("Earth heliocentric:\t"); printVect(v_e);
    printf("      geocentric:\t");
    VECT3_DIFF(v_g, v_e, v_e); printVect(v_g);

    printBody(BODY_MERCURY, "Mercury", T, v_e, e);
    printBody(BODY_VENUS, "Venus", T, v_e, e);
    printBody(BODY_MARS, "Mars", T, v_e, e);
    printBody(BODY_JUPITER, "Jupiter", T, v_e, e);
    printBody(BODY_SATURN, "Saturn", T, v_e, e);
    printBody(BODY_URANUS, "Uranus", T, v_e, e);

    vsop87_clearDataCache();

    return 0;
}

