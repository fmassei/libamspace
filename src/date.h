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
#ifndef H_DATE_H
#define H_DATE_H

#include "exports.h"
#include "precision.h"
#include <time.h>

/* like a struct tm, but see struct description */
typedef struct dtm_s dtm_st;

/* struct tm is fine, but tm_sec is an int, while we need more precision:
 * thus we will use the struct tm members except for the seconds */
struct dtm_s {
    struct tm tm;
    DT sec;
};

/* dtm-julian conversions */
LIB_PUBLIC DT amsp_date_tm2julian(const dtm_st *dtm);
LIB_PUBLIC void amsp_date_YD2tm(dtm_st *out, int year, DT days);
LIB_PUBLIC void amsp_date_julian2tm(dtm_st *out, DT j);
/* get greenwich time */
LIB_PUBLIC DT amsp_date_gstime(DT j);

#endif /* H_DATE_H */
