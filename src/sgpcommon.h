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
#ifndef H_SGPCOMMON_H
#define H_SGPCOMMON_H

#include "precision.h"

typedef struct sat_s sat_st;
typedef struct orbitalelems_s orbitalelems_st;

/* classical orbital elements */
struct orbitalelems_s {
    DT mm;          /* (n(_o)) mean motion (at epoch) */
    DT e;           /* (e(_o)) mean eccentricity (at epoch) */
    DT i;           /* (i(_o)) mean inclination (at epoch) */
    DT ma;          /* (M(_o)) mean anomaly (at epoch) */
    DT argp;        /* (w(_o)) mean argument of perigee (at epoch) */
    DT node;        /* (omega(_o)) mean longitude of ascending node (at epoch)*/
};

/* error in computation */
typedef enum saterr_e {
    SATERR_NOERR=0, /* everything went fine */
    SATERR_MM,      /* inconsistent mm has been calculated */
    SATERR_E,       /* inconsistent e " " " */
    SATERR_PL,      /* inconsistent pl " " " */
    SATERR_DECAY,   /* decay condition */
} saterr_et;


#endif /* H_SGPCOMMON_H */
