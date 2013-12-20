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
#ifndef H_SGPUTILS_H
#define H_SGPUTILS_H

#include "exports.h"
#include "precision.h"
#include "sat.h"
#include "gravconst.h"

LIB_LOCAL void utils_unkozai(sat_st *sat, const gravconst_st *g);
LIB_LOCAL void utils_getDensityConstants(const sat_st *sat,
                                    const gravconst_st *g, DT *ss, DT* qms4);

#endif /* H_SGPUTILS_H */
