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
#ifndef H_TLE_H
#define H_TLE_H

#include "exports.h"
#include <stdio.h>

typedef struct tle_s tle_st;

struct tle_s {
    char name[25];
    int satNo;
    char class;
    int designatorYY;
    int designatorLaunchNo;
    char designatorPiece[4];
    int epochYY;
    double epoch;
    double dtmm;    /* first time derivative of mean motion / 2 */
    double dt2mm;   /* second time derivative of mean motion / 6 */
    double bstar;   /* BSTAR drag term */
    int elSet;
    double inclinationDeg;
    double rightAscDeg;
    double ecc;
    double argOfPerigeeDeg;
    double meanAnomalyDeg;
    double meanMotionRpD;
    int revs;
    tle_st *next;
};

LIB_PUBLIC void tle_destroy(tle_st **tle);
/* read one or more TLEs from file
 * readAll: 1 = will read the whole file until EOF and return the TLE list
 *          0 = read the first TLE and return
 * compatMode: 1 = compability mode, skip comments and don't expect header
 *             0 = expect header+line1+line2, no comments
 */
LIB_PUBLIC tle_st *tle_read(FILE *in, int readAll, int compatMode);

#endif /* H_TLE_H */
