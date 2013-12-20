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
#include <functions.h>
#include <legendre.h>
#include <interpolate.h>
#include <integrate.h>
#include <linear.h>
#include <leastSquare.h>

#include <stdio.h>
#include <stdlib.h>

/* sine function */
static DT fSin(DT x) { return SIN(x); }
/* Legendre polinomial of order 10 */
static DT fLeg10(DT x) { return amsp_legendrePoly(10, x); }
/* sin(x)*cos(2x)^(3/2) */
static DT fFnc1(DT x) { return POW(COS(2.*x),1.5)*SIN(x); }

/* general vector checker */
#define ISNEAR(_X, _Y) (FABS((_X)-(_Y))<1e-5)
static int matchRes(DT *v1, DT *v2, int n)
{
    int i;
    for (i=0; i<n; ++i)
        if (!ISNEAR(v1[i], v2[i]))
            return 0;
    return 1;
}
#undef ISNEAR

/* **************************************************** functions */
static int test_functions(void)
{
    DT testRes[] = { 0., PI, PI*2, PI*3 };
    DT out[4];
    int nSol;
    if ((nSol = amsp_func1r_findZeros(fSin, 0, 4*PI, 500, out, 4))!=4)
        return EXIT_FAILURE;
    if (!matchRes(out, testRes, 4))
        return EXIT_FAILURE;
    return 0;
}

/* **************************************************** legendre */
static int test_legendre(void)
{
    DT testRes[] = { -0.973907, -0.865063, -0.67941, -0.433395, -0.148874,
                      0.148874,  0.433395,  0.67941,  0.865063,  0.973907 };
    DT out[10];
    int nSol;
    if ((nSol = amsp_func1r_findZeros(fLeg10, -1, 1, 5000, out, 10))!=10)
        return EXIT_FAILURE;
    if (!matchRes(out, testRes, 10))
        return EXIT_FAILURE;
    return 0;
}

/* **************************************************** interpolations */
static int test_interpolate1(void)
{
    DT testRes[] = { .23, -0.049275, -0.0075555 };
    DT out[3];
    amsp_interpolate_getBesselCoeff(.73, out);
    if (!matchRes(out, testRes, 3))
        return EXIT_FAILURE;
    return 0;
}
static int test_interpolate2(void)
{
    DT testRes[] = {.121445};
    DT toInterpolate[] = { .381300, .285603, .190092,
                           .096327, .008268, -.067725 },
       at = .073/.1;
    DT out[1];
    out[0] = amsp_interpolate_bessel(toInterpolate, 2, at);
    if (!matchRes(out, testRes, 1))
        return EXIT_FAILURE;
    return 0;
}
static int test_interpolate3(void)
{
    DT testRes[] = { .756221, .800624 };
    DT toInterpolate[] = { 0., .5, .86603, 1. },
       at = 21./30.;
    DT out[2];
    out[0] = amsp_interpolate_linear(toInterpolate[1], toInterpolate[2], at);
    out[1] = amsp_interpolate_bessel(toInterpolate, 1, at);
    if (!matchRes(out, testRes, 2))
        return EXIT_FAILURE;
    return 0;
}
static int test_interpolate(void)
{
    if (    test_interpolate1() ||
            test_interpolate2() ||
            test_interpolate3() )
        return EXIT_FAILURE;
    return 0;
}

/* **************************************************** integrations */
static int test_integrate(void)
{
    DT testRes[] = { 1.00228, 1.00228, 0.108769 };
    DT x[11], y[11], qCoeff[3];
    DT out[3];
    amsp_func1r_tabulate(fSin, 0., PI/2., x, y, 3);
    out[0] = amsp_integrate_simpsonRule(x[0], y[0], x[1], y[1], x[2], y[2]);
    amsp_interpolate_3pointsQuadratic(x[0], y[0], x[1], y[1], x[2], y[2],
                                      qCoeff);
    out[1] = amsp_integrate_quadratic(qCoeff, PI/4, PI/4);
    amsp_func1r_tabulate(fFnc1, 0., PI/4., x, y, 11);
    out[2] = amsp_integrate_simpson(x, y, 11);
    if (!matchRes(out, testRes, 3))
        return EXIT_FAILURE;
    return 0;
}

/* **************************************************** linear */
static int test_linear_solve(void)
{
    DT testRes[] = { 2., 7., 6., 4., 3. };
    DT out[5];
    int nOut;
    DT e0[] = { 9, -9, 8, -6, 4, -9 },
       e1[] = { 5, -1, 6, 1, 5, 58 },
       e2[] = { 2, 4, -5, -6, 7, -1 },
       e3[] = { 2, 3, -8, -5, -2, -49 },
       e4[] = { 8, -5, 7, 1, 5, 42 },
       *sys[] = { e0, e1, e2, e3, e4 };
    if (amsp_linear_solve(sys, 5, 6, out, &nOut)!=0)
        return EXIT_FAILURE;
    if (!matchRes(out, testRes, 5))
        return EXIT_FAILURE;
    return 0;
}
static int test_linear_solveResidue(void)
{
    DT testRes[] = { 2.473933, 5.396512, 3.722554 };
    DT out[5];
    int nOut;
    DT e0[] = { 7, -6, 8, 15 },
       e1[] = { 3, 5, -2, 27 },
       e2[] = { 2, -2, 7, 20 },
       e3[] = { 4, 2, -5, 2 },
       e4[] = { 9, -8, 7, 5 },
       *sys[] = { e0, e1, e2, e3, e4 };
    if (amsp_linear_solveResidue(sys, 5, 4, out, &nOut)!=0)
        return EXIT_FAILURE;
    if (!matchRes(out, testRes, 3))
        return EXIT_FAILURE;
    return 0;
}
static int test_linear(void)
{
    if (test_linear_solve()!=0 || test_linear_solveResidue()!=0)
        return EXIT_FAILURE;
    return 0;
}

/* **************************************************** least squares */
static int test_leastSquareLine(void)
{
    DT testRes0[] = { .55, .9 },
       testRes1[] = { 1.5493, -.950704 };
    DT x[] = { 1, 2, 3, 4, 5 },
       y[] = { 1, 2.5, 2.75, 3, 3.5 };
    DT out[2];
    if (amsp_leastSquare_coeff(x, y, 5, &out[0], &out[1])!=0)
        return EXIT_FAILURE;
    if (!matchRes(out, testRes0, 2))
        return EXIT_FAILURE;
    if (amsp_leastSquare_coeff(y, x, 5, &out[0], &out[1])!=0)
        return EXIT_FAILURE;
    if (!matchRes(out, testRes1, 2))
        return EXIT_FAILURE;
    return 0;
}
static int test_leastSquareQuad(void)
{
    DT testRes[] = { -0.002247, 3.774837, -961.336719 };
    DT x[] = { 395.1, 448.1, 517.7, 583.3, 790.2 },
       y[] = { 171.0, 289.0, 399.0, 464.0, 620.0 };
    DT out[3];
    if (amsp_leastSquare_quad(x, y, 5, out)!=0)
        return EXIT_FAILURE;
    printf("%f %f %f\n", out[0], out[1], out[2]);
    if (!matchRes(out, testRes, 3))
        return EXIT_FAILURE;
    return 0;
}
static int test_leastsquare(void)
{
    if (test_leastSquareLine()!=0 || test_leastSquareQuad()!=0)
        return EXIT_FAILURE;
    return 0;
}

/* **************************************************** main */
int main(const int argc, const char *argv[])
{
    if (argc!=2)
        return EXIT_FAILURE;
    if (*argv[1]=='f') {
        return test_functions();
    } else if (*argv[1]=='l') {
        return test_legendre();
    } else if (*argv[1]=='i') {
        return test_interpolate();
    } else if (*argv[1]=='n') {
        return test_integrate();
    } else if (*argv[1]=='e') {
        return test_linear();
    } else if (*argv[1]=='s') {
        return test_leastsquare();
    } else {
        fprintf(stderr, "Unknown test %c\n", *argv[1]);
        return EXIT_FAILURE;
    }
}

