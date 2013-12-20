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
#ifndef H_EXPORTS_H
#define H_EXPORTS_H

#if defined _WIN32 || defined __CYGWIN__
#   ifdef BUILDING_LIBAMSPACE
#       ifdef __GNUC__
#           define LIB_PUBLIC __attribute__ ((dllexport))
#       else
#           define LIB_PUBLIC __declspec(dllexport)
#       endif
#   else
#       ifdef __GNUC__
#           define LIB_PUBLIC __attribute__ ((dllimport))
#       else
#           define LIB_PUBLIC __declspec(dllimport)
#       endif
#   endif
#   define LIB_LOCAL
#else
#   if __GNUC__ >= 4
#       define LIB_PUBLIC __attribute__ ((visibility ("default")))
#       define LIB_LOCAL __attribute__ ((visibility ("hidden")))
#   else
#       define LIB_PUBLIC
#       define LIB_LOCAL
#   endif
#endif

#endif /* H_EXPORTS_H */
