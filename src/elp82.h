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
#ifndef H_ELP82_H
#define H_ELP82_H

#include "exports.h"
#include "precision.h"

typedef struct elp82data_s elp82data_st;

/* free all the ELP82 data */
LIB_PUBLIC void amsp_elp82_destroyData(elp82data_st **data);

/* load ELP82 data */
LIB_PUBLIC elp82data_st *amsp_elp82_getData(const char *datadir);
/* get moon position at date */
LIB_PUBLIC void amsp_elp82_getPos(const elp82data_st *data, double jdate,
                                  double out[3]);

#endif /* H_ELP82_H */
