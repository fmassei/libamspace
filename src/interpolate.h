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
#ifndef H_INTERPOLATE_H
#define H_INTERPOLATE_H

#include "exports.h"
#include "precision.h"

LIB_PUBLIC DT amsp_interpolate_linear(DT v1, DT v2, DT rt);
LIB_PUBLIC void amsp_interpolate_getBesselCoeff(DT v, DT *out);
LIB_PUBLIC DT amsp_interpolate_bessel(const DT *val, int from, DT rt);

#endif /* H_INTERPOLATE_H */
