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
#ifndef H_VSOP87FILE_H
#define H_VSOP87FILE_H

#include "exports.h"
#include "precision.h"

typedef struct vsop87_hdr_s vsop87_hdr_st;
typedef struct vsop87_rcd_s vsop87_rcd_st;
typedef struct vsop87_filepart_s vsop87_filepart_st;

struct vsop87_hdr_s {
    int versionCode;    /* 0 = VSOP87 (main), 1=A, 2=B, ..., 5=E */
    char name[8];       /* EARTH, MARS, etc. */
    int coordIndex;
    int alpha;
    int nTerms;
};
struct vsop87_rcd_s {
    int versionCode;
    int bodyCode;
    int coordIndex;
    int alpha;
    int n;              /* rank in serie */
    int a[12];          /* coefficients of mean longitudes */
    DT S;
    DT K;
    DT A;
    DT B;
    DT C;
};
struct vsop87_filepart_s {
    vsop87_hdr_st hdr;
    vsop87_rcd_st *rcds;
    vsop87_filepart_st *next;
};

typedef enum vsop87_body_e {
    BODY_MERCURY=1,
    BODY_VENUS,
    BODY_EARTH,
    BODY_MARS,
    BODY_JUPITER,
    BODY_SATURN,
    BODY_URANUS,
    BODY_NEPTUNE,
    BODY_EARTHMOON,
} vsop87_body_et;

LIB_LOCAL void vsop87_filepart_destroy(vsop87_filepart_st **fp);
LIB_LOCAL vsop87_filepart_st *vsop87_filepart_read(const char *datadir,
                                                vsop87_body_et body, char code);

#endif /* H_VSOP87FILE_H */
