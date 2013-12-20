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
#ifndef H_INTEGRATE_H
#define H_INTEGRATE_H

#include "exports.h"
#include "precision.h"

LIB_PUBLIC void amsp_interpolate_3pointsQuadratic(DT x0, DT y0, DT x1, DT y1,
                                                  DT x2, DT y2, DT *qCoeff);

LIB_PUBLIC DT amsp_integrate_quadratic(const DT *qCoeff, DT x, DT dx);
LIB_PUBLIC DT amsp_integrate_simpsonRule(DT x0, DT y0, DT x1, DT y1,
                                         DT x2, DT y2);
LIB_PUBLIC DT amsp_integrate_simpson(const DT *x, const DT *y, int n);

#endif /* H_INTEGRATE_H */
