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
#include "linear.h"
#include <stdlib.h>

static void reverse(DT *arr, int nArr)
{
    int i, j;
    DT sw;
    for (i=0, j=nArr-1; i<j; ++i, --j) {
        sw = arr[i];
        arr[i] = arr[j];
        arr[j] = sw;
    }
}
static int solveHelper(DT **eqs, int nEq, int nCf, DT *sol, int *nSol)
{
    DT **step;
    DT nt;
    int i, ret = 0;
    if (nEq<=1)
        return 0;
    if ((step = calloc((nEq-1), sizeof(*step)))==NULL)
        return -1;
    for (i=1; i<nEq; ++i) {
        DT c1 = eqs[0][0],
           c2 = eqs[i][0];
        int j;
        if ((step[i-1] = calloc((nCf-1), sizeof(**step)))==NULL) {
            ret = -1;
            goto freeandexit;
        }
        for (j=1; j<nCf; ++j)
            step[i-1][j-1] = c2*eqs[0][j]-c1*eqs[i][j];
    }
    if ((ret = solveHelper(step, nEq-1, nCf-1, sol, nSol))!=0)
        goto freeandexit;
    nt = step[0][nCf-2];
    for (i=1; i<nCf-2; ++i)
        nt -= step[0][i]*sol[*nSol-i];
    sol[*nSol] = nt/step[0][0];
    *nSol += 1;
freeandexit:
    for (i=0; i<nEq-1; ++i)
        free(step[i]);
    free(step);
    return ret;
}

int amsp_linear_solve(DT **eqs, int nEq, int nCf, DT *sol, int *nSol)
{
    DT nt;
    int i;
    *nSol = 0;
    if (solveHelper(eqs, nEq, nCf, sol, nSol)!=0)
        return -1;
    nt = eqs[0][nCf-1];
    for (i=1; i<nCf-1; ++i)
        nt -= eqs[0][i]*sol[*nSol-i];
    sol[*nSol] = nt/eqs[0][0];
    *nSol += 1;
    reverse(sol, *nSol);
    return 0;
}


int amsp_linear_solveResidue(DT **eqs, int nEq, int nCf, DT *sol, int *nSol)
{
    DT **equiv;
    int i, ret;
    if ((equiv = calloc(nCf-1, sizeof(*equiv)))==NULL)
        return -1;
    for (i=0; i<nCf-1; ++i) {
        int j;
        if ((equiv[i] = calloc(nCf, sizeof(**equiv)))==NULL) {
            ret = -1;
            goto freeandexit;
        }
        for (j=0; j<nCf; ++j) {
            DT coeff = 0.;
            int k;
            for (k=0; k<nEq; ++k)
                coeff += eqs[k][i]*eqs[k][j];
            equiv[i][j] = coeff;
        }
    }
    ret = amsp_linear_solve(equiv, nCf-1, nCf, sol, nSol);
freeandexit:
    for (i=0; i<nCf-1; ++i)
        free(equiv[i]);
    free(equiv);
    return ret;
}

