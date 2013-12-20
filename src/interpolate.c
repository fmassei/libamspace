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
#include "interpolate.h"

DT amsp_interpolate_linear(DT v1, DT v2, DT rt)
{
    return rt*v2+(1-rt)*v1;
}

void amsp_interpolate_getBesselCoeff(DT v, DT *out)
{
    out[0] = v-0.5;
    out[1] = v*(v-1)/4.;
    out[2] = v*(0.5+v*(-1.5+v))/6.;
    out[3] = v*(2.+v*(-1+v*(-1+v)))/48.;
    out[4] = v*(-1.+v*(2.5+v*v*(-2.5+v)))/120;
}

DT amsp_interpolate_bessel(const DT *vals, int from, DT rt)
{
    return -(1/6.)*rt*((2-rt*(3-rt))*vals[from-1] +
            (1-rt)*vals[from+2]) +
            (1/2.)*((2+rt*(-1+rt*(-2+rt)))*vals[from] +
                    rt*(2+rt*(1-rt))*vals[from+1]);
}

