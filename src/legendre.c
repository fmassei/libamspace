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
#include "legendre.h"

DT amsp_legendrePoly(int n, DT v)
{
    DT p=0, pn1, pn2, k;
    int i;
    if (n==0) return 1.;
    if (n==1) return v;
    pn1 = 1; pn2 = 0;
    for (i=2; i<=n+1; ++i) {
        k = (DT)(i-1.);
        p = ((2*k-1)/k)*v*pn1-((k-1)/k)*pn2;
        pn2 = pn1;
        pn1 = p;
    }
    return p;
}

