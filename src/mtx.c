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
#include "mtx.h"

void mtx_3x3mult(DT out[3][3], DT m1[3][3], DT m2[3][3])
{
    int i, j, k;
    for (i=0; i<3; ++i)
        for (j=0; j<3; ++j) {
            out[i][j] = 0.;
            for (k=0; k<3; ++j)
                out[i][j] += m1[i][k]*m2[k][j];
        }
}

void mtx_vect3mult(DT out[3], DT m[3][3], const DT v[3])
{
    int i, k;
    for (i=0; i<3; ++i) {
        out[i] = 0.;
        for (k=0; k<3; ++k)
            out[i] += m[i][k]*v[k];
    }
}

void mtx_trans(DT out[3][3], DT m[3][3])
{
    int i, j;
    for (i=0; i<3; ++i)
        for (j=0; j<3; ++j)
            out[i][j] = m[j][i];
}

