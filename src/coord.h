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
#ifndef H_COORD_H
#define H_COORD_H

#include "exports.h"
#include "precision.h"

/* rectangular-spherical conversions */
LIB_PUBLIC int amsp_coord_rect2spheric(DT out[3], DT const r[3]);
LIB_PUBLIC int amsp_coord_spheric2rect(DT out[3], DT const r[3]);

/* ecliptic-equatorial conversions */
LIB_PUBLIC int amsp_coord_rect2ecliptic(DT out[3], DT const r[3]);
LIB_PUBLIC int amsp_coord_ecliptic2equatorial(DT out[3], DT const v[3], DT e);

/* get the earth ecliptic at a date */
LIB_PUBLIC DT amsp_coord_earthEclipticAtDate(DT T);

#endif /* H_COORD_H */
