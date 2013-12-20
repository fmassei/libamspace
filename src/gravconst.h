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
#ifndef H_GRAVCONST_H
#define H_GRAVCONST_H

#include "exports.h"
#include "precision.h"

typedef struct gravconst_s gravconst_st;

typedef enum gravconsttype_e {
    WGS72OLD,           /* we use this just to test old test-cases */
    WGS72,
    WGS84
} gravconsttype_et;

struct gravconst_s {
    /* datum-dependent */
    DT aE;              /* = 1. */
    DT rE;              /* Earth radius (Km) */
    DT ke;              /* (ke) sqrt(G/M)(er/min)^(3/2) */
    DT j2;              /* 2nd gravitation zonal harmonic of Earth */
    DT j3;              /* 3rd gravitation zonal harmonic of Earth */
    DT j4;              /* 4th gravitation zonal harmonic of Earth */
    DT prE;             /* Earth polar radius (Km) */
    DT f;               /* flattening */
    /* sun / moon constants */
    DT cL, cS;          /* .25*(mx/mx+me)*Nx (mx=body mass, me=earth mass, Nx=
                         * apparent mm of body (sun=1, moon=me*1/81.53) */
    DT iLe;             /* moon inclination on ecliptic */
    DT nL, nS;          /* (mean apparent motion?) */
    DT eL, eS;          /* (apparent eccentricity?) */
    /* pre-calculated */
    DT k2;              /* (1/2)j2(a_e^2) */
    DT k4;              /* (-3/8)j4(a_e^4) */
    DT a30cof;          /* -j3/k2*(a_e^3) */
    DT s, qms4;         /* parameters for the SGP4/SGP8 density function */
};

/* TODO: change name! */
LIB_PUBLIC void getgravconst(gravconsttype_et type, gravconst_st *out);

#endif /* H_GRAVCONST_H */
