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
#include <stdio.h>
#include <stdlib.h>
#include <sgp4.h>
#include <coord.h>
#include <projections.h>

int main(void)
{
    tle_st *tle;
    gravconst_st g;
    DT start, stop, delta;

    getgravconst(WGS84, &g);
    start = 0.; stop = 1440; delta = 10;
    while ((tle = tle_read(stdin, 0, 0))!=NULL) {
        sat_st sat;
        DT sph[3], mer[2], utm[4];
        satFromTLE(&sat, tle);
        sgp4_init(&sat, &g);
        sat.t = start;
        while (sat.t<=stop) {
            int ret;
            if ((ret = sgp4(&sat, &g))!=0) {
                fprintf(stderr, "error %d\n", sat.error);
                goto nextsat;
            }
            amsp_coord_rect2spheric(sph, sat.pos);
            amsp_projections_sph2latlon(sph, sph);
            /*sph[1] += (360./(24.*60))*sat.t;
            if (sph[1]>180.) sph[1]-=360.;*/
            amsp_projections_latlon2mercator(mer, sph, 2048, 1588);
            amsp_projections_latlon2utm(utm, sph, &g);
            printf("%f %f > %f %f > %f %f %f\n",
                sph[0], sph[1], mer[0], mer[1], utm[0], utm[1], utm[2]);
            /*printf("%f %f\n",
                mer[0], mer[1]);*/
            sat.t += delta;
        }
nextsat:
        tle_destroy(&tle);
    }
    return 0;
}

