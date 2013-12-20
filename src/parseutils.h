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
#ifndef H_PARSEUTILS_H
#define H_PARSEUTILS_H

#include "exports.h"
#include <stdlib.h>

LIB_LOCAL void parseutils_buf2str(char *out, const char *buf, size_t len);
LIB_LOCAL int parseutils_buf2int(int *out, const char *buf, size_t len);
LIB_LOCAL int parseutils_buf2double(double *out, const char *buf, size_t len);
LIB_LOCAL int parseutils_buf2doubleDPA(double *out, const char *buf, size_t len);

#endif /* H_PARSEUTILS_H */
