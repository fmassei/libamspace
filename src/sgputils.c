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
#include "sgputils.h"
#include "precision.h"

/* from the sat structure, calculate and assign the original mean motion (omm),
 * the semimajor axis (sma) and the perigee */
void utils_unkozai(sat_st *sat, const gravconst_st *g)
{
    DT a1, d1, a0, d0;
    a1 = POW(g->ke/sat->o.mm, 2./3.);
    d1 = (3.*g->k2*(3.*P2(COS(sat->o.i))-1.)) /
         (2.*P2(a1)*POW(1.-P2(sat->o.e), 1.5));
    a0 = a1*(1.-(1./3.)*d1-P2(d1)-(134./81.)*P3(d1));
    d0 = (3.*g->k2*(3.*P2(COS(sat->o.i))-1.)) /
         (2.*P2(a0)*POW(1.-P2(sat->o.e), 1.5));
    sat->omm = sat->o.mm/(1.+d0);
    sat->sma = a0/(1.-d0);
    sat->perigee = (sat->sma*(1-sat->o.e)-g->aE)*g->rE;
    sat->isDeep = ((2.*PI/sat->omm)>=225.);
}

void utils_getDensityConstants(const sat_st *sat, const gravconst_st *g,
                               DT *ss, DT* qms4)
{
    if (sat->perigee<156.) {
        if (sat->perigee<98.) {
            *ss = 20.;
        } else {
            //*ss = sat->sma*(1.-sat->o.e)-g->s;
            *ss = sat->perigee-78.;
        }
        //*qms4 = P4(POW(g->qms4, .25)+g->s-*ss);
        *qms4 = POW(((120.-*ss)*g->aE/g->rE),4);
        *ss = *ss/g->rE+g->aE;
    } else {
        *ss = g->s;
        *qms4 = g->qms4;
    }
}

