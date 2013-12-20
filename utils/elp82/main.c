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
#include <elp82.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void printUsage(void)
{
    fprintf(stderr, "Usage: elp82 <datadir> <juliandate>\n"
                    "\tif juliandate==\"test\" run the internal test\n");
}
static inline void printVect(DT v[3])
{
    printf("%.9f %.9f %.9f\n", v[0], v[1], v[2]);
}

int main(const int argc, const char *argv[])
{
    const char *datadir;
    DT jdate;
    elp82data_st *data;
    DT out[3];

    if (argc!=3) {
        printUsage();
        return EXIT_FAILURE;
    }
    datadir = argv[1];
    if ((data = amsp_elp82_getData(datadir))==NULL)
        return EXIT_FAILURE;
    if (strcmp(argv[2], "test")!=0) {
        jdate = strtod(argv[2], NULL);
        amsp_elp82_getPos(data, jdate, out);
        printVect(out);
    } else {
        DT res[5][4] = {
            { 2469000.5, -361602.98535627,  44996.99510257, -30696.65315716 },
            { 2449000.5, -363132.34247954,  35863.65378335, -33196.00409288 },
            { 2429000.5, -371577.58161204,  75271.14315490, -32227.94618056 },
            { 2409000.5, -373896.15893047, 127406.79128901, -30037.79225087 },
            { 2389000.5, -346331.77361283, 206365.40364249, -28502.11731834 } };
        int i;
        for (i=0; i<5; ++i) {
            amsp_elp82_getPos(data, res[i][0], out);
            if (    FABS(out[0]-res[i][1])>1e-8 ||
                    FABS(out[1]-res[i][2])>1e-8 ||
                    FABS(out[2]-res[i][3])>1e-8 ) {
                fprintf(stderr, "Failed on test #%d.\n", i);
                return EXIT_FAILURE;
            }
        }
    }
    amsp_elp82_destroyData(&data);
    return 0;
}

