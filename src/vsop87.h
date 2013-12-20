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
#ifndef H_VSOP87_H
#define H_VSOP87_H

#include "exports.h"
#include "vsop87file.h"
#include "precision.h"

typedef enum vsop87type_e {
    VSOP87_HELIO_ELL_J2000,
    VSOP87_HELIO_RECT_J2000,
    VSOP87_HELIO_SPH_J2000,
    VSOP87_HELIO_RECT_DATE,
    VSOP87_HELIO_SPH_DATE,
    VSOP87_BARY_RECT_J2000,
} vsop87type_et;

LIB_PUBLIC int vsop87_getCoords(DT out[6], const char *datadir,
                            vsop87type_et type, vsop87_body_et body, DT date);
LIB_PUBLIC void vsop87_clearDataCache(void);

#endif /* H_VSOP87_H */
