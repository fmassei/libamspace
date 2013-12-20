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
#include "leastSquare.h"
#include "linear.h"
#include <stdlib.h>

int amsp_leastSquare_coeff(DT *x, DT *y, int n, DT *m, DT *q)
{
    DT **eqs, sol[2];
    int i, ret = 0, nSol;
    if ((eqs = calloc(n, sizeof(*eqs)))==NULL)
        return -1;
    for (i=0; i<n; ++i) {
        if ((eqs[i] = calloc(3, sizeof(**eqs)))==NULL) {
            ret = -1;
            goto freeandexit;
        }
        eqs[i][0] = x[i];
        eqs[i][1] = 1;
        eqs[i][2] = y[i];
    }
    ret = amsp_linear_solveResidue(eqs, n, 3, sol, &nSol);
    *m = sol[0];
    *q = sol[1];
freeandexit:
    for (i=0; i<n; ++i)
        free(eqs[i]);
    free(eqs);
    return ret;
}

int amsp_leastSquare_quad(DT *x, DT *y, int n, DT *qCoef)
{
    DT **eqs;
    int i, ret = 0, nQcoef;
    if ((eqs = calloc(n, sizeof(*eqs)))==NULL)
        return -1;
    for (i=0; i<n; ++i) {
        if ((eqs[i] = calloc(4, sizeof(**eqs)))==NULL) {
            ret = -1;
            goto freeandexit;
        }
        eqs[i][0] = P2(x[i]);
        eqs[i][1] = x[i];
        eqs[i][2] = 1;
        eqs[i][3] = y[i];
    }
    ret = amsp_linear_solveResidue(eqs, n, 4, qCoef, &nQcoef);
freeandexit:
    for (i=0; i<n; ++i)
        free(eqs[i]);
    free(eqs);
    return ret;
}

