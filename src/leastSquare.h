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
#ifndef H_LEASTSQUARE_H
#define H_LEASTSQUARE_H

#include "exports.h"
#include "precision.h"

/* least squares linear */
LIB_PUBLIC int amsp_leastSquare_coeff(DT *x, DT *y, int n, DT *m, DT *q);
/* least squares quadratic */
LIB_PUBLIC int amsp_leastSquare_quad(DT *x, DT *y, int n, DT *qCoef);

#endif /* H_LEASTSQUARE_H */
