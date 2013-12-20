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
#include "angles.h"
#include <stdio.h>

void amsp_angle_deg2sex(sexangle_st *sex, DT deg)
{
    sex->deg = (int)FLOOR(deg);
    sex->min = (int)(FMOD(FABS(deg*60.), 60.));
    sex->sec = FMOD(FABS(deg*3600.), 60.);
}
void amsp_angle_sex2deg(DT *deg, const sexangle_st *sex)
{
    DT rs = sex->deg*3600+sex->min*60+sex->sec;
    *deg = rs/3600.;
}

void amsp_angle_sexprint(const sexangle_st *sex)
{
    printf("%d°%d'"PF, sex->deg, sex->min, sex->sec);
}

