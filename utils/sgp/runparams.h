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
#ifndef H_RUNPARAMS_H
#define H_RUNPARAMS_H

#include <gravconst.h>

typedef struct runparams_s runparams_st;
typedef enum runalg_e {
    RUNALG_SGP4,
    RUNALG_SGP8,
} runalg_et;

struct runparams_s {
    FILE *in, *out;                     /* IO streams */
    gravconst_st gravconst;             /* gravitational constants */
    double startmfe, stopmfe, deltamin; /* time interval */
    runalg_et runalg;                   /* which algorithm to run */
    int TLEcompabilityMode;             /* if reading TLE in compability mode */
    int testModeOutput;                 /* test mode output */
};

#endif /* H_RUNPARAMS_H */
