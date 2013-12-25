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
#ifndef H_NUTATIONS_H
#define H_NUTATIONS_H

#include "exports.h"
#include "precision.h"

typedef struct nutdata_s nutdata_st;

struct nutdata_s {
    int iar[107][5];
    DT rar[107][4];
};

LIB_PUBLIC nutdata_st *amsp_nutation_getData(const char *fname);
LIB_PUBLIC void amsp_nutations_freeData(nutdata_st **nut);

#endif /* H_NUTATIONS_H */
