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
#ifndef H_PRECISION_H
#define H_PRECISION_H

#ifndef LONGDOUBLES
#   include <math.h>
#   define DT       double

#   define PI       M_PI
#   define COS      cos
#   define SIN      sin
#   define TAN      tan
#   define ASIN     asin
#   define ACOS     acos
#   define ATAN2    atan2
#   define POW      pow
#   define SQRT     sqrt
#   define FMOD     fmod
#   define FABS     fabs
#   define STRTOD   strtod
#   define FLOOR    floor
#   define LOG      log

#   define PF       "%f"

#else

#   define __USE_GNU
#   include <math.h>
#   define DT       long double

#   define PI       M_PIl
#   define COS      cosl
#   define SIN      sinl
#   define TAN      tanl
#   define ASIN     asinl
#   define ACOS     acosl
#   define ATAN2    atan2l
#   define POW      powl
#   define SQRT     sqrtl
#   define FMOD     fmodl
#   define FABS     fabsl
#   define STRTOD   strtold
#   define FLOOR    floorl
#   define LOG      logl

#   define PF       "%Lf"
#endif

#define ZERO        1e-12
#define ISZERO(_X)  (FABS(_X)<ZERO)

/* powers, from ^2 to ^8 */
#define P2(_X)      ((_X)*(_X))
#define P3(_X)      (P2(_X)*(_X))
#define P4(_X)      (P3(_X)*(_X))
#define P5(_X)      (P4(_X)*(_X))
#define P6(_X)      (P5(_X)*(_X))
#define P7(_X)      (P6(_X)*(_X))
#define P8(_X)      (P7(_X)*(_X))

/* euclidean distance */
#define DIST2(_X, _Y)       (SQRT(P2(_X)+P2(_Y)))
#define DIST3(_X, _Y, _Z)   (SQRT(P2(_X)+P2(_Y)+P2(_Z)))

/* vector operations */
#define VECT3_SUM(_OUT, _V0, _V1)   do {(_OUT)[0] = (_V0)[0]+(_V1)[0];\
                                        (_OUT)[1] = (_V0)[1]+(_V1)[1];\
                                        (_OUT)[2] = (_V0)[2]+(_V1)[2];} while(0)
#define VECT3_DIFF(_OUT, _V0, _V1)  do {(_OUT)[0] = (_V0)[0]-(_V1)[0];\
                                        (_OUT)[1] = (_V0)[1]-(_V1)[1];\
                                        (_OUT)[2] = (_V0)[2]-(_V1)[2];} while(0)
#define VECT3_SCALAR(_OUT, _S, _V0) do {(_OUT)[0] = (_S)*(_V0)[0];\
                                        (_OUT)[1] = (_S)*(_V0)[1];\
                                        (_OUT)[2] = (_S)*(_V0)[2];} while(0)

#endif /* H_PRECISION_H */
