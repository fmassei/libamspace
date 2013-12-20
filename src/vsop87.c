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
#include "vsop87.h"
#include <stdlib.h>

static vsop87_filepart_st *s_byBody[6][9] = {{0}};

static char codeFromType(vsop87type_et type)
{
    switch (type) {
    case VSOP87_HELIO_RECT_J2000: return 'A';
    case VSOP87_HELIO_SPH_J2000: return 'B';
    case VSOP87_HELIO_RECT_DATE: return 'C';
    case VSOP87_HELIO_SPH_DATE: return 'D';
    case VSOP87_BARY_RECT_J2000: return 'E';
    case VSOP87_HELIO_ELL_J2000: default: return ' ';
    }
}

static vsop87_filepart_st *vsop87_getData(const char *datadir,
                                      vsop87_body_et body, vsop87type_et type)
{
    if (s_byBody[type][body-1]==NULL)
        s_byBody[type][body-1] = vsop87_filepart_read(datadir, body,
                                                      codeFromType(type));
    return s_byBody[type][body-1];
}

int vsop87_getCoords(DT out[6], const char *datadir,
                            vsop87type_et type, vsop87_body_et body, DT date)
{
    vsop87_filepart_st *fp, *p;
    DT Ta[6];
    unsigned i;
    if ((fp = vsop87_getData(datadir, body, type))==NULL)
        return -1;
    out[0] = out[1] = out[2] = out[3] = out[4] = out[5] = 0.;
    date = (date-2451545.)/365250.;
    Ta[0] = 1.; Ta[1] = date;
    for (i=2; i<6; ++i)
        Ta[i] = Ta[i-1]*Ta[1];
    for (p=fp; p!=NULL; p=p->next) {
        int i;
        for (i=0; i<p->hdr.nTerms; ++i) {
            DT u, cosu;
            vsop87_rcd_st *r = &p->rcds[i];
            u = r->B+r->C*Ta[1];
            cosu = COS(u);
            out[r->coordIndex-1] += Ta[r->alpha]*r->A*cosu;
            if (type!=VSOP87_HELIO_ELL_J2000) {
                DT pt = (r->alpha>0) ? (Ta[r->alpha-1]*r->alpha*r->A*cosu) : 0;
                out[r->coordIndex+2] += pt-Ta[r->alpha]*r->A*r->C*SIN(u);
            }
        }
    }
    if (type!=VSOP87_HELIO_ELL_J2000) {
        out[3] /= 365250.;
        out[4] /= 365250.;
        out[5] /= 365250.;
    }
    if (type==VSOP87_HELIO_ELL_J2000) {
        out[1] = FMOD(out[1], 2.*PI);
        if (out[1]<0.) out[1] += 2.*PI;
    } else if (type==VSOP87_HELIO_SPH_J2000 || type==VSOP87_HELIO_SPH_DATE) {
        out[0] = FMOD(out[0], 2.*PI);
        if (out[0]<0.) out[0] += 2.*PI;
    }
    return 0;
}

void vsop87_clearDataCache(void)
{
    unsigned i, j;
    for (i=0; i<6; ++i)
        for (j=0; j<9; ++j)
            if (s_byBody[i][j])
                vsop87_filepart_destroy(&s_byBody[i][j]);
}

