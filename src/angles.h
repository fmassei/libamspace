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
#ifndef H_ANGLES_H
#define H_ANGLES_H

#include "exports.h"
#include "precision.h"

#define CGR (PI/180.)       /* Convert deGrees to Radians */
#define CRG (180./PI)       /* Convert Radians to deGrees */

/* sexagesimal angles */
typedef struct sexangle_s sexangle_st;

struct sexangle_s {
    int deg, min;
    DT sec;
};

/* convert between degrees and sexagesimal */
LIB_PUBLIC void amsp_angle_deg2sex(sexangle_st *sex, DT deg);
LIB_PUBLIC void amsp_angle_sex2deg(DT *deg, const sexangle_st *sex);

/* print out a sexagesimal angle */
LIB_PUBLIC void amsp_angle_sexprint(const sexangle_st *sex);

#endif /* H_ANGLES_H */
