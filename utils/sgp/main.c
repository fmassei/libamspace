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
#include <getopt.h>
#include <libamspace.h>
#include <sgp4.h>
#include <sgp8.h>
#include "runparams.h"

static void printUsage(void)
{
    fprintf(stderr,
"Usage: sgp [options]\n"
"\t-h, --help                           this message\n"
"\t-i, --infile <fname>                 input file (default: stdin)\n"
"\t-o, --outfile <fname>                output file (default: stdout)\n"
"\t-S, --start <min>                    start time (default: 0)\n"
"\t-s, --stop <min>                     stop time (default: 1440)\n"
"\t-d, --delta <min>                    delta time (default: 360)\n"
"\t-D, --datum <WGS72OLD|WGS72|WGS84>   datum (default:WGS84)\n"
"\t-m, --method <4|8>                   select SGP4 or SGP8 (default: SGP8)\n"
"\t-c, --compat                         read the TLE in compability mode\n"
"\t-t, --test                           test mode output\n"
"\t-v, --version                        print version and exit\n"
"\n");
}

static int parseOptions(const int argc, char * const *argv, runparams_st *rp)
{
    struct option long_options[] = {
        { "help", 0, 0, 'h' },
        { "infile", 1, 0, 'i' },
        { "outfile", 1, 0, 'o' },
        { "start", 1, 0, 'S' },
        { "stop", 1, 0, 's' },
        { "delta", 1, 0, 'd' },
        { "datum", 1, 0, 'D' },
        { "method", 1, 0, 'm' },
        { "compat", 0, 0, 'c' },
        { "test", 0, 0, 't' },
        { "version", 0, 0, 'v' },
        { 0, 0, 0, 0 } };
    rp->in = stdin;
    rp->out = stdout;
    rp->startmfe = 0.;
    rp->stopmfe = 1440.;
    rp->deltamin = 360.;
    rp->runalg = RUNALG_SGP8;
    getgravconst(WGS84, &rp->gravconst);
    rp->TLEcompabilityMode = 0;
    rp->testModeOutput = 0;
    while(1) {
        int c, option_index = 0;
        c = getopt_long(argc, argv, "hi:o:S:s:d:D:m:ctv",
                        long_options, &option_index);
        if (c==-1)
            break;
        switch (c) {
        case 'h':
            printUsage();
            exit(EXIT_SUCCESS);
            break;
        case 'i':
            if ((rp->in = fopen(optarg, "rb"))==NULL) {
                perror("fopen");
                return -1;
            }
            break;
        case 'o':
            if ((rp->out = fopen(optarg, "wb"))==NULL) {
                perror("fopen");
                return -1;
            }
            break;
        case 'S':
            rp->startmfe = STRTOD(optarg, NULL);
            break;
        case 's':
            rp->stopmfe = STRTOD(optarg, NULL);
            break;
        case 'd':
            rp->deltamin = STRTOD(optarg, NULL);
            break;
        case 'D':
            if (!strcmp(optarg, "WGS72OLD"))
                getgravconst(WGS72OLD, &rp->gravconst);
            else if (!strcmp(optarg, "WGS72"))
                getgravconst(WGS72, &rp->gravconst);
            else if (!strcmp(optarg, "WGS84"))
                getgravconst(WGS84, &rp->gravconst);
            else {
                fprintf(stderr, "Unknown datum %s\n", optarg);
                return -1;
            }
            break;
        case 'm': {
            int i;
            i = atoi(optarg);
            if (i==4) rp->runalg = RUNALG_SGP4;
            else if (i==8) rp->runalg = RUNALG_SGP8;
            else {
                fprintf(stderr, "Unknown method %d\n", i);
                return -1;
            }
            }
            break;
        case 'c':
            rp->TLEcompabilityMode = 1;
            break;
        case 't':
            rp->TLEcompabilityMode = 1;
            rp->startmfe = -1440.;
            rp->stopmfe = 1440.;
            rp->deltamin = 10.;
            rp->runalg = RUNALG_SGP4;
            getgravconst(WGS72, &rp->gravconst);
            rp->testModeOutput = 1;
            break;
        case 'v':
            fprintf(stderr, "libamspace ver %d.%d.%s\n",
                LIBAMSPACE_VERSION_MAJOR, LIBAMSPACE_VERSION_MINOR,
                LIBAMSPACE_VERSION_REV);
            exit(EXIT_SUCCESS);
            break;
        case '?':
        case ':':
        default:
            fprintf(stderr, "Invalid options.\n\n");
            printUsage();
            return -1;
        }
    }
    if (optind<argc) {
        fprintf(stderr, "Invalid extra parameters\n\n");
        printUsage();
        return -1;
    }
    return 0;
}

int main(const int argc, char * const *argv)
{
    tle_st *tle;
    runparams_st rp;
    int ret;
    if (parseOptions(argc, argv, &rp)!=0)
        return EXIT_FAILURE;
    while ((tle = tle_read(rp.in, 0, rp.TLEcompabilityMode))!=NULL) {
        sat_st sat;
        satFromTLE(&sat, tle);
        if (rp.testModeOutput)
            fprintf(rp.out, "%d xx\n", tle->satNo);
        if (rp.runalg==RUNALG_SGP4) sgp4_init(&sat, &rp.gravconst);
        else sgp8_init(&sat, &rp.gravconst);
        
        if (rp.testModeOutput) {
            sat.t = 0.;
            if (rp.runalg==RUNALG_SGP4) ret = sgp4(&sat, &rp.gravconst);
            else ret = sgp8(&sat, &rp.gravconst);
            if (ret!=0) {
                fprintf(stderr, "%d error %d\n", tle->satNo, sat.error);
                goto nextsat;
            }
            fprintf(rp.out, " %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f\n",
                    sat.t, sat.pos[0], sat.pos[1], sat.pos[2],
                    sat.vel[0], sat.vel[1], sat.vel[2]);
        }
        sat.t = rp.startmfe;
        while (sat.t<=rp.stopmfe) {
            dtm_st dtm;
            DT jd;
            if (rp.runalg==RUNALG_SGP4) ret = sgp4(&sat, &rp.gravconst);
            else ret = sgp8(&sat, &rp.gravconst);
            if (ret!=0) {
                fprintf(stderr, "%d error %d\n", tle->satNo, sat.error);
                goto nextsat;
            }
            jd = sat.epochJD+sat.t/1440.;
            amsp_date_julian2tm(&dtm, jd);
            fprintf(rp.out, " %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f %04d/%02d/%02d %02d:%02d:%9.6f\n",
                sat.t, sat.pos[0], sat.pos[1], sat.pos[2],
                       sat.vel[0], sat.vel[1], sat.vel[2],
                       dtm.tm.tm_year+1900, dtm.tm.tm_mon+1, dtm.tm.tm_mday,
                       dtm.tm.tm_hour, dtm.tm.tm_min, dtm.sec);
            sat.t += rp.deltamin;
            //printf("* i=%f e=%f a=%f n=%f\n", sat.c.i*CRG, sat.c.e, sat.c.argp*CRG, sat.c.node*CRG);
            //ln(tan(PI/4+lat/2.))
        }
        /*fprintf(stderr, "%d end: %s", tle->satNo, sat.isDeep?"deep":"nearth");
        if (sat.isDeep)
            fprintf(stderr, " res %s", sat.deep.res==RESONANCE_NONE?"no":
                                   (sat.deep.res==RESONANCE_ONEDAY?"one":"half"));
        fprintf(stderr, " %lf\n", sat.bstar);*/
nextsat:
        tle_destroy(&tle);
    }
    return 0;
}

