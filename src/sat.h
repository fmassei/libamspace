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
#ifndef H_SAT_H
#define H_SAT_H

#include "exports.h"
#include "precision.h"
#include "tle.h"
#include "sgpcommon.h"
#include "deep.h"

typedef struct nearth_s nearth_st;
/* near-earth variables */
struct nearth_s {
    int isVeryNear; /* perigee less than 220Km */
    DT eta;
    DT D[3];        /* secular coefficients */
    DT dw, dm;      /* secular effects of atmospheric drag and gravitation 
                     * for argp and ma*/
    DT dmo;         /* secular gravity on ma */
};

struct sat_s {
    /* orbital quantities */
    orbitalelems_st o;  /* orbital quantities at epoch */
    /* derivatives */
    DT mmdt;            /* (ndo) first derivative of mean motion at epoch */
    DT mmd2t;           /* (nddo) second derivative of mean motion at epoch */
    /* drag coefficient */
    DT bstar;           /* SGP4 drag coefficient */
    /* epoch */
    DT epochJD;         /* epoch in julian days */

    /* general, precalculated quantities */
    DT omm;             /* (n''_o) original mean motion */
    DT sma;             /* (a''_o) semimajor axis */
    DT perigee;         /* perigee */
    DT C[5];            /* drag constants */
    DT Tcof[4];         /* time coefficients */
    DT xcof, ycof;      /* long periodic coefficients */
    DT ncof;            /* secular gravity on node coefficient */
    DT mdt, adt, ndt;   /* SGP8 only: ma, argp and node partial increment */
    DT gamma, ed, pp, qq, ovgpp; /* SGP8 only: ?? TODO */

    int isDeep;         /* 1 if in "deep space" */
    nearth_st nearth;   /* near Earth data */
    dspace_st deep;     /* deep space data */

    /* current status (updated each iteration) */
    saterr_et error;    /* error, set on premature exit */
    DT t;               /* time (minutes) */
    orbitalelems_st c;  /* "current" orbital quantities */
    orbitalelems_st c_d;/* current derivatives */
    DT d_omm;           /* SGP8 only, n_0' */
    DT pos[3];          /* (r) position (km) */
    DT vel[3];          /* (rd) velocity (km/s) */
};

LIB_PUBLIC void satFromTLE(sat_st *sat, const tle_st *tle);

#endif /* H_SAT_H */
