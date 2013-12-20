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
#ifndef H_FUNCTIONS_H
#define H_FUNCTIONS_H

#include "exports.h"
#include "precision.h"

typedef DT(*func1r_ft)(DT);

LIB_PUBLIC void amsp_func1r_tabulate(func1r_ft f, DT from, DT to,
                                     DT *x, DT *y, int Nout);
LIB_PUBLIC int amsp_func1r_bisect(func1r_ft f, DT a, DT b, int maxIter,
                                     DT *out);
LIB_PUBLIC int amsp_func1r_findZeros(func1r_ft f, DT a, DT b, int nStep,
                                     DT *out, int Nout);

#endif /* H_FUNCTIONS_H */
