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
#ifndef H_MTX_H
#define H_MTX_H

#include "exports.h"
#include "precision.h"

LIB_LOCAL void mtx_3x3mult(DT out[3][3], DT m1[3][3], DT m2[3][3]);
LIB_LOCAL void mtx_vect3mult(DT out[3], DT m[3][3], const DT v[3]);
LIB_LOCAL void mtx_trans(DT out[3][3], DT m[3][3]);

#endif /* H_MTX_H */
