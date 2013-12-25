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
#include "nutations.h"
#include <stdio.h>
#include <stdlib.h>

static nutdata_st *s_nutdata = NULL;

nutdata_st *amsp_nutation_getData(const char *fname)
{
    nutdata_st *ret;
    FILE *in;
    int i, j;
    if (s_nutdata!=NULL)
        return s_nutdata;
    if ((in = fopen(fname, "rb"))==NULL)
        return NULL;
    if ((ret = malloc(sizeof(*ret)))==NULL) {
        fclose(in);
        return NULL;
    }
    for (i=0; i<107; ++i) {
        if (fscanf(in, "%d %d %d %d %d %lf %lf %lf %lf %d\n",
                    &ret->iar[i][0], &ret->iar[i][1], &ret->iar[i][2],
                    &ret->iar[i][3], &ret->iar[i][4],
                    &ret->rar[i][0], &ret->rar[i][1], &ret->rar[i][2],
                    &ret->rar[i][3], &j)==EOF) {
            amsp_nutations_freeData(&ret);
            fclose(in);
            return NULL;
        }
        for (j=0; j<4; ++j)
            ret->rar[i][j] *= .0001/3600.;
    }
    fclose(in);
    s_nutdata = ret;
    return ret;
}

void amsp_nutations_freeData(nutdata_st **nut)
{
    if (nut==NULL || *nut==NULL) return;
    free(*nut);
    *nut = s_nutdata = NULL;
}

