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
#include "integrate.h"

void amsp_interpolate_3pointsQuadratic(DT x0, DT y0, DT x1, DT y1, DT x2, DT y2,
                                       DT *qCoeff)
{
    DT dx = FABS(x1-x0);
    qCoeff[0] = y1 - ((x1*(y2-y0))/(2*dx)) +
                ((P2(x1)*(y2-2*y1+y0))/(2*P2(dx)));
    qCoeff[1] = ((y2-y0)/(2*dx)) - ((x1*(y2-2*y1+y0))/P2(dx));
    qCoeff[2] = (y2-2*y1+y0)/(2*P2(dx));
    return;
    x2=x2; /* unused */
}

DT amsp_integrate_quadratic(const DT *qCoeff, DT x, DT dx)
{
    return 2*dx*(qCoeff[0]+qCoeff[1]*x+qCoeff[2]*P2(x)+(1/3.)*qCoeff[2]*P2(dx));
}

DT amsp_integrate_simpsonRule(DT x0, DT y0, DT x1, DT y1, DT x2, DT y2)
{
    DT dx = FABS(x1-x0);
    return ((1/3.)*(y0+4*y1+y2)*dx);
    x2=x2; /* unused */
}

DT amsp_integrate_simpson(const DT *x, const DT *y, int n)
{
    DT ret = 0.;
    int i;
    for (i=0; i<n-1; i+=2) {
        if (i+2>=n)
            return ret;
        ret += amsp_integrate_simpsonRule(x[i], y[i], x[i+1], y[i+1],
                                          x[i+2], y[i+2]);
    }
    return ret;
}

