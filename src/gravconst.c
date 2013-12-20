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
#include "gravconst.h"
#include "angles.h"
#include <math.h>

void getgravconst(gravconsttype_et type, gravconst_st *out)
{
    double mu;
    switch(type) {
    case WGS72OLD:
        out->rE = 6378.135;
        out->ke = 0.0743669161;
        out->j2 = 0.001082616;
        out->j3 = -0.00000253881;
        out->j4 = -0.00000165597;
        out->f = 298.26;
        out->prE = 6356.7505;
        break;
    case WGS72:
        mu = 398600.8;
        out->rE = 6378.135;
        out->ke = 60./SQRT(P3(out->rE)/mu);
        out->j2 = 0.001082616;
        out->j3 = -0.00000253881;
        out->j4 = -0.00000165597;
        out->f = 298.26;
        out->prE = 6356.7505;
        break;
    case WGS84:
        mu = 398600.5;
        out->rE = 6378.137;
        out->ke = 60./SQRT(P3(out->rE)/mu);
        out->j2 = 0.00108262998905;
        out->j3 = -0.00000253215306;
        out->j4 = -0.00000161098761;
        out->f = 298.257223563;
        out->prE = 6356.7523142;
        break;
    }

    out->aE = 1.;

    out->cL = 4.7968065e-7;
    out->cS = 2.98647969e-6;
    out->iLe = 5.145396*CGR;
    out->nL = 1.5835218e-4;
    out->nS = 1.19459e-5;
    out->eL = .05490;
    out->eS = .01675;

    out->k2 = .5*out->j2;
    out->k4 = (-3./8.)*out->j4;
    out->a30cof = -out->j3/out->k2*P3(out->aE);
    out->s = 78./out->rE + 1.;
    out->qms4 = P4((120.-78.)/out->rE);
}

