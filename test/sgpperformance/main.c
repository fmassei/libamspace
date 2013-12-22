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
#include <sys/time.h>
#include <pthread.h>

static gravconst_st g;

static double elapsedTime(void) { // returns 0 seconds first time called
    static struct timeval t0;
    struct timeval tv;
    gettimeofday(&tv, 0);
    if (!t0.tv_sec)
        t0 = tv;
    return tv.tv_sec - t0.tv_sec + (tv.tv_usec - t0.tv_usec) / 1000000.;
}
static int test_nonThreaded(tle_st *tle)
{
    tle_st *p;
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
            sat.t += 1;
        }
    }
    return 0;
}

static void *threadP(tle_st *p)
{
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
        sat.t += 1;
    }
    return NULL;
}
static int test_threaded(tle_st *tle)
{
    tle_st *p;
    pthread_t th[15000];
    int i=0;
    for (p=tle; p!=NULL; p=p->next, ++i) {
        pthread_create(&th[i], NULL, threadP, p);
    }
    for (; i>=0; --i)
        pthread_join(th[i], NULL);
    return 0;
}

pthread_mutex_t mt;
static void *threadP2(tle_st *tle)
{
    tle_st *p;
    for (p=tle; p!=NULL; p=p->next) {
        sat_st sat;
        pthread_mutex_lock(&mt);
        if (p->dtmm<-1111111.) {
            pthread_mutex_unlock(&mt);
            continue;
        }
        p->dtmm = -200000000.;
        pthread_mutex_unlock(&mt);
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
            sat.t += 1;
        }
    }
    return NULL;
}
static void test_threaded2(tle_st *tle)
{
    pthread_t th[16];
    int i;
    for (i=0; i<16; ++i)
        pthread_create(&th[i], NULL, threadP2, tle);
    for (i=0; i<16; ++i)
        pthread_join(th[i], NULL);
    return 0;
}
int main(const int argc, const char *argv[])
{
    tle_st *tle;
    FILE *in;
    double et;
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
    test_nonThreaded(tle);
    et = elapsedTime()-et;
    fprintf(stderr, "elapsed time: %f\n", et);

    et = elapsedTime();
    test_threaded(tle);
    et = elapsedTime()-et;
    fprintf(stderr, "elapsed time: %f\n", et);

    et = elapsedTime();
    test_threaded2(tle);
    et = elapsedTime()-et;
    fprintf(stderr, "elapsed time: %f\n", et);

    tle_destroy(&tle);
    return 0;
}

