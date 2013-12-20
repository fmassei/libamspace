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
#include <libamspace.h>
#include <linear.h>
#include <stdio.h>
#include <stdlib.h>

int main(void)
{
    DT **eqs = NULL;
    int nEq, nCf, i, j;
    char mod;
    DT *sol = NULL;
    int nSol = 0;

    printf("Linear equation solver.\n");
    printf("[N]ormal or [R]esidue: "); fflush(stdout);
    scanf("%c", &mod);
    if (mod!='N' && mod!='R') {
        fprintf(stderr, "Unknown method %c\n", mod);
        return EXIT_FAILURE;
    }
    printf("Number of equations: "); fflush(stdout);
    scanf("%d", &nEq);
    printf("Number of variables: "); fflush(stdout);
    scanf("%d", &nCf);
    if ((eqs = calloc(nEq, sizeof(*eqs)))==NULL) {
        perror("calloc");
        return EXIT_FAILURE;
    }
    for (i=0; i<nEq; ++i) {
        if ((eqs[i] = calloc(nCf+1, sizeof(**eqs)))==NULL) {
            perror("calloc");
            return EXIT_FAILURE;
        }
        for (j=0; j<nCf+1; ++j)
            scanf("%lf", &eqs[i][j]);
    }

    if ((sol = calloc(nCf, sizeof(*sol)))==NULL) {
        perror("calloc");
        return EXIT_FAILURE;
    }
    if (mod=='N' && nEq==nCf) {
        printf("Solving system:\n");
        if (amsp_linear_solve(eqs, nEq, nCf+1, sol, &nSol)!=0) {
            fprintf(stderr, "error solving\n");
            return EXIT_FAILURE;
        }
        for (i=0; i<nSol; ++i)
            printf("x%d = %f\n", i, sol[i]);
    } else if (nEq>nCf) {
        printf("Solving with residue:\n");
        if (amsp_linear_solveResidue(eqs, nEq, nCf+1, sol, &nSol)!=0) {
            fprintf(stderr, "error solving\n");
            return EXIT_FAILURE;
        }
        for (i=0; i<nSol; ++i)
            printf("x%d = %f\n", i, sol[i]);
    } else {
        fprintf(stderr, "Too few equations.\n");
    }
    for (i=0; i<nEq; ++i)
        free(eqs[i]);
    free(eqs);
    free(sol);
    return 0;
}

