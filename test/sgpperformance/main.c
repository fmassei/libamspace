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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libamspace.h>
#include <sgp4.h>

#include <time.h>

double elapsedTime(void) { // returns 0 seconds first time called
    static struct timeval t0;
    struct timeval tv;
    gettimeofday(&tv, 0);
    if (!t0.tv_sec)
        t0 = tv;
    return tv.tv_sec - t0.tv_sec + (tv.tv_usec - t0.tv_usec) / 1000000.;
}
int main(const int argc, const char *argv[])
{
    tle_st *tle, *p;
    FILE *in;
    double et;
    gravconst_st g;
    if (argc!=2) {
        fprintf(stderr, "Usage: sgpperformance <TLEfile>\n");
        return EXIT_FAILURE;
    }
    if ((in = fopen(argv[1], "rb"))==NULL) {
        perror("fopen");
        return EXIT_FAILURE;
    }
    if ((tle = tle_read(in, 1, 0))==NULL) {
        fprintf(stderr, "Error reading TLEs\n");
        return EXIT_FAILURE;
    }
    fclose(in);
    getgravconst(WGS84, &g);
    et = elapsedTime();
    for (p=tle; p!=NULL; p=p->next) {
        sat_st sat;
        satFromTLE(&sat, p);
        sgp4_init(&sat, &g);
        sat.t = 0;
        while (sat.t<=1440) {
            if (sgp4(&sat, &g)!=0) {
                fprintf(stderr, "Error %d on %d\n", p->satNo, sat.error);
                break;
            }
            /*printf("%d %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f\n",
                p->satNo, sat.t,
                          sat.pos[0], sat.pos[1], sat.pos[2],
                          sat.vel[0], sat.vel[1], sat.vel[2]);*/
            sat.t += 10;
        }
    }
    et = elapsedTime()-et;
    fprintf(stderr, "elapsed time: %f\n", et);
    tle_destroy(&tle);
    return 0;
}

