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
#include "functions.h"

void amsp_func1r_tabulate(func1r_ft f, DT from, DT to, DT *x, DT *y, int Nout)
{
    DT step, i;
    int k;
    step = (to-from)/(DT)(Nout-1);
    for (k=0, i=from; k<Nout; ++k, i+=step)
        y[k] = f(x[k]=i);
}

static inline int sign(DT x)
{
    return (x>0)-(x<0);
}
#define TOL 1e-7    /* TODO: use a parameter instead? */
int amsp_func1r_bisect(func1r_ft f, DT a, DT b, int maxIter, DT *out)
{
    DT c;
    int i;
    for (i=0; i<maxIter; ++i) {
        c = (a+b)/2.;
        if (FABS(f(c))<TOL || ((b-a)/2.)<TOL) {
            *out = c;
            return 0;
        }
        if (sign(f(c))==sign(f(a))) a = c;
        else b = c;
    }
    return -1;
}

int amsp_func1r_findZeros(func1r_ft f, DT a, DT b, int nStep,
                                                            DT *out, int Nout)
{
    DT xStep = (b-a)/(DT)nStep,
       x;
    int sgn, i;
    for (sgn=sign(f(a)), x=a, i=0; x<b && i<Nout; x+=xStep) {
        if (sign(f(x))!=sgn) {
            DT zero;
            if (amsp_func1r_bisect(f, a, x, nStep, &zero)!=0)
                return -1;
            out[i++] = zero;
            a = x;
            sgn = sign(f(a));
        }
    }
    return i;
}

