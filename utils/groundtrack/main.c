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
    start = 0.; stop = 1440; delta = .1;
    while ((tle = tle_read(stdin, 0, 0))!=NULL) {
        sat_st sat;
        DT sph[3], mer[2], utm[4];
        satFromTLE(&sat, tle);
        sgp4_init(&sat, &g);

        {
            dtm_st dtm;
            struct tm *tm;
            time_t t;
            DT jdate;
            while (1) {
            t = time(NULL);
            tm = gmtime(&t);
            dtm.tm = *tm;
            dtm.sec = tm->tm_sec;
            jdate = amsp_date_tm2julian(&dtm);
            printf("%f\n", jdate);
            //jdate -= amsp_date_gstime(jdate);
            //printf("%f\n", jdate);
            printf("%f\n", sat.epochJD);
            sat.t = (jdate-sat.epochJD)*1440;
            printf("%f\n", sat.t);
            sgp4(&sat, &g);
            amsp_coord_rect2spheric(sph, sat.pos);
            sph[1] -= amsp_date_gstime(jdate);
            amsp_coord_ecliptic2equatorial(sph, sph, amsp_coord_earthEclipticAtDate(jdate));
            amsp_projections_sph2latlon(sph, sph);
            printf("%f %f\n", sph[1], sph[0]);
            /*amsp_coord_rect2spheric(sph, sat.pos);
            amsp_coord_ecliptic2equatorial(sph, sph, amsp_coord_earthEclipticAtDate(jdate));
            amsp_projections_sph2latlon(sph, sph);
            printf("%f %f\n", sph[1], sph[0]);*/
            sleep(1);
            }
        }
        continue;

        sat.t = start;
        while (sat.t<=stop) {
            int ret;
            if ((ret = sgp4(&sat, &g))!=0) {
                fprintf(stderr, "error %d\n", sat.error);
                goto nextsat;
            }
            amsp_coord_rect2spheric(sph, sat.pos);
            /*printf("%f %f %f\n",
                sat.pos[0], sat.pos[1], sat.pos[2]);*/
            amsp_projections_sph2latlon(sph, sph);
            //printf("%f %f\n",
            //    sph[1], sph[0]);
            sph[1] += (360./(24.*60))*sat.t;
            if (sph[1]>180.) sph[1]-=360.;
            printf("%f %f\n",
                sph[1], sph[0]);
            amsp_projections_latlon2mercator(mer, sph, 2048, 1588);
            amsp_projections_latlon2utm(utm, sph, &g);
            /*printf("%f %f > %f %f > %f %f %f > %f\n",
                sph[1], sph[0], mer[0], mer[1], utm[0], utm[1], utm[2],
                (360./(24.*60))*sat.t);*/
            /*printf("%f %f\n",
                mer[0], mer[1]);*/
            //printf("%f %f\n", sph[1], utm[1]);
            sat.t += delta;
        }
nextsat:
        tle_destroy(&tle);
    }
    return 0;
}

