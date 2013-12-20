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
#include <string.h>
#include "sat.h"
#include "date.h"
#include "angles.h"

#define XPDOTP  (1440./(2*PI))

void satFromTLE(sat_st *sat, const tle_st *tle)
{
    dtm_st dtm;
    memset(sat, 0, sizeof(*sat));
    amsp_date_YD2tm(&dtm, (tle->epochYY<57) ? tle->epochYY+2000
                                            : tle->epochYY+1900, tle->epoch);
    sat->epochJD = amsp_date_tm2julian(&dtm);
    sat->mmdt = tle->dtmm/(XPDOTP*1440.);
    sat->mmd2t = tle->dt2mm/(XPDOTP*P2(1440.));
    sat->bstar = tle->bstar;
    sat->o.i = tle->inclinationDeg*CGR;
    sat->o.node = tle->rightAscDeg*CGR;
    sat->o.e = tle->ecc;
    sat->o.argp = tle->argOfPerigeeDeg*CGR;
    sat->o.ma = tle->meanAnomalyDeg*CGR;
    sat->o.mm = tle->meanMotionRpD/XPDOTP;
}

