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
#include "coord.h"
#include "angles.h"

/* r = sqrt(x^2+y^2+z^2)
 * cos(theta) = z/r
 * tan(ro) = y/x
*/
int amsp_coord_rect2spheric(DT out[3], DT const r[3])
{
    /*out[0] = ATAN2(r[1], r[0]);
    out[1] = ATAN2(r[2], DIST2(r[0], r[1]));*/
    out[0] = (PI/2.)-ACOS(r[1]/DIST3(r[0], r[1], r[2]));
    out[1] = ATAN2(r[0],r[2]);
    out[2] = DIST3(r[0], r[1], r[2]);
    return 0;
}
/* x = rsin(theta)cos(ro)
 * y = rsin(theta)sin(ro)
 * z = rcos(theta) */
int amsp_coord_spheric2rect(DT out[3], DT const r[3])
{
    out[0] = r[2]*SIN(out[0])*COS(out[1]);
    out[0] = r[2]*SIN(out[0])*SIN(out[1]);
    out[0] = r[2]*COS(out[0]);
    return 0;
}

/* tan(lon) = Y/X
 * cos(lat) = Z/sqrt(X^2+y^2)
 * dist = sqrt(X^2+Y^2+Z^2) */
int amsp_coord_rect2ecliptic(DT out[3], DT const r[3])
{
    return amsp_coord_rect2spheric(out, r);
}

/* ecliptic: lon=g lat=b
 * equatorial: rightasc=a decl=d
 * tana = (singcose-tanbsine)/cosg
 * sind = sinbcose+cosbsinesing */
int amsp_coord_ecliptic2equatorial(DT out[3], DT const v[3], DT e)
{
    DT g=v[0], b=v[1];
    out[0] = ATAN2(SIN(g)*COS(e)-TAN(b)*SIN(e), COS(g));
    out[1] = ASIN(SIN(b)*COS(e)+COS(b)*SIN(e)*SIN(g));
    out[2] = v[2];
    return 0;
}

DT amsp_coord_earthEclipticAtDate(DT T)
{
    sexangle_st trms[] = { {23, 26, 21.406},
                           { 0,  0, 46.836769},
                           { 0,  0,  0.0001831},
                           { 0,  0,  0.00200340},
                           { 0,  0,  0.576e-6},
                           { 0,  0,  4.34e-8} };
    int signs[] = { 1, -1, -1, 1, -1, -1 };
    DT ret = 0;
    unsigned i;
    T = (T-2451545.)/365250.;
    for (i=0; i<sizeof(trms)/sizeof(trms[0]); ++i) {
        DT fct;
        amsp_angle_sex2deg(&fct, &trms[i]);
        ret += signs[i]*fct*POW(T,i);
    }
    return ret*CGR;
}

