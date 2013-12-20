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
#include <string.h>
#include <time.h>

static const char *s_datadir;

#define ISNEAR(_X, _Y) (FABS((_X)-(_Y))<1e-4)
static int matchRes(DT *v1, DT *v2, int n)
{
    int i;
    for (i=0; i<n; ++i)
        if (!ISNEAR(v1[i], v2[i]))
            return 0;
    return 1;
}

static void printVect(DT *v, int n)
{
    int i;
    for (i=0; i<n; ++i) printf("%f ", v[i]);
    printf("\n");
}

static void printUsage(void)
{
    fprintf(stderr, "Usage: vsop87test <datadir> <testfile>\n");
}

static int parseHdr(const char *buf, vsop87type_et* type, vsop87_body_et *body,
                    DT *T)
{
    switch(buf[0]) {
    case '_': *type = VSOP87_HELIO_ELL_J2000; break;
    case 'A': *type = VSOP87_HELIO_RECT_J2000; break;
    case 'B': *type = VSOP87_HELIO_SPH_J2000; break;
    case 'C': *type = VSOP87_HELIO_RECT_DATE; break;
    case 'D': *type = VSOP87_HELIO_SPH_DATE; break;
    case 'E': *type = VSOP87_BARY_RECT_J2000; break;
    default: return -1;
    }
    if (!memcmp(buf+2, "MER", 3)) *body = BODY_MERCURY;
    else if (!memcmp(buf+2, "VEN", 3)) *body = BODY_VENUS;
    else if (!memcmp(buf+2, "MAR", 3)) *body = BODY_MARS;
    else if (!memcmp(buf+2, "JUP", 3)) *body = BODY_JUPITER;
    else if (!memcmp(buf+2, "SAT", 3)) *body = BODY_SATURN;
    else if (!memcmp(buf+2, "URA", 3)) *body = BODY_URANUS;
    else if (!memcmp(buf+2, "NEP", 3)) *body = BODY_NEPTUNE;
    else if (!memcmp(buf+2, "EMB", 3)) *body = BODY_EARTHMOON;
    else if (!memcmp(buf+2, "EAR", 3)) *body = BODY_EARTH;
    else if (!memcmp(buf+2, "SUN", 3)) *body = BODY_EARTHMOON;
    *T = strtod(buf+6, NULL);
    return 0;
}
int main(const int argc, const char *argv[])
{
    FILE *in;
    char buf[200];
    int ret = 0, line = 0;
    if (argc!=3) {
        printUsage();
        return EXIT_FAILURE;
    }
    s_datadir = argv[1];
    if ((in = fopen(argv[2], "rb"))==NULL) {
        perror("fopen");
        return EXIT_FAILURE;
    }
    while (fgets(buf, sizeof(buf), in)!=NULL) {
        DT our[6], tst[6];
        vsop87type_et type;
        vsop87_body_et body;
        DT T;
        ++line;
        if (buf[strlen(buf)-1]!='\n' && buf[strlen(buf)-1]!='\r') {
            ret = EXIT_FAILURE;
            fprintf(stderr, "Error reading line %d\n", line);
            goto freeandexit;
        }
        if (parseHdr(buf, &type, &body, &T)!=0) {
            ret = EXIT_FAILURE;
            fprintf(stderr, "Error reading line %d\n", line);
            goto freeandexit;
        }
        if (vsop87_getCoords(our, s_datadir, type, body, T)!=0) {
            ret = EXIT_FAILURE;
            fprintf(stderr, "Error in vsop87_getCoords, line %d\n", line);
            goto freeandexit;
        }
        if (    fgets(buf, sizeof(buf), in)==NULL ||
                (buf[strlen(buf)-1]!='\n' && buf[strlen(buf)-1]!='\r') ) {
            ret = EXIT_FAILURE;
            fprintf(stderr, "Error reading line %d\n", line);
            goto freeandexit;
        }
        ++line;
        sscanf(buf, "%lf,%lf,%lf,%lf,%lf,%lf",
                &tst[0], &tst[1], &tst[2], &tst[3], &tst[4], &tst[5]);
        switch (type) {
        case VSOP87_HELIO_ELL_J2000:
            if (    !ISNEAR(our[0], tst[0]) ||      /* a */
                    !ISNEAR(our[1], tst[3]) ||      /* l */
                    !ISNEAR(our[2], tst[1]) ||      /* k */
                    !ISNEAR(our[3], tst[4]) ||      /* h */
                    !ISNEAR(our[4], tst[2]) ||      /* q */
                    !ISNEAR(our[5], tst[5]) ) {     /* p */
                ret = EXIT_FAILURE;
                fprintf(stderr, "Test failed at line %d\n", line);
                printVect(our, 6);
                printf("%lf %lf %lf %lf %lf %lf\n",
                    tst[0], tst[3], tst[1], tst[4], tst[2], tst[5]);
                goto freeandexit;
            }
            break;
        case VSOP87_HELIO_RECT_J2000:
        case VSOP87_HELIO_RECT_DATE:
        case VSOP87_BARY_RECT_J2000:
            if (!matchRes(our, tst, 6)) {
                ret = EXIT_FAILURE;
                fprintf(stderr, "Test failed at line %d\n", line);
                printVect(our, 6);
                printVect(tst, 6);
                goto freeandexit;
            }
            break;
        case VSOP87_HELIO_SPH_J2000:
        case VSOP87_HELIO_SPH_DATE:
            if (!matchRes(our, tst, 6)) {
                ret = EXIT_FAILURE;
                fprintf(stderr, "Test failed at line %d\n", line);
                printVect(our, 6);
                printVect(tst, 6);
                goto freeandexit;
            }
            break;
        }
    }
freeandexit:
    fclose(in);
    vsop87_clearDataCache();
    return ret;
}

