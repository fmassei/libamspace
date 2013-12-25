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
#include "fk5.h"
#include "nutations.h"
#include "mtx.h"
#include "angles.h"

static void precess(DT out[3][3], DT jcent)
{
    DT cv = PI/(180.*3600.),
       zeta, theta, z,
       czeta, szeta, ctheta, stheta, cz, sz;
    zeta = ((.017998*jcent+.30188)*jcent+2306.2181)*jcent;
    theta = ((-.041833*jcent-.42665)*jcent+2004.3109)*jcent;
    z = ((.018203*jcent+1.09468)*jcent+2306.2181)*jcent;
    czeta = COS(zeta*cv); szeta = SIN(zeta*cv);
    ctheta = COS(theta*cv); stheta = SIN(theta*cv);
    cz = COS(z*cv); sz = SIN(z*cv);
    out[0][0] = czeta*ctheta*cz - szeta*sz;
    out[0][1] = czeta*ctheta*sz + szeta*cz;
    out[0][2] = czeta*stheta;
    out[1][0] = -szeta*ctheta*cz - czeta*sz;
    out[1][1] = -szeta*ctheta*sz + czeta*cz;
    out[1][2] = -szeta*stheta;
    out[2][0] = -stheta*cz;
    out[2][1] = -stheta*sz;
    out[2][2] = ctheta;
}

static void truemean(DT out[3][3], DT jcent, const nutdata_st *nut)
{
    DT nt[3][3], st[3][3];
    DT meps, dpsi, deps, teps, eqe;
    DT cpsi, spsi, ceps, seps, cteps, steps;
    DT l, l1, f, d, o;
    int i;
    
    meps = ((.001813*jcent-.00059)*jcent-46.8150)*jcent+84381.448;
    meps = FMOD(meps/3600.,360.)*CGR;
    
    l = (((.064*jcent+31.310)*jcent+1717915922.6330)*jcent)/3600.+134.96298139;
    l1 = (((-.012*jcent-.577)*jcent+129596581.2240)*jcent)/3600.+357.52772333;
    f = (((.011*jcent-13.257)*jcent+1739527263.1370)*jcent)/3600.+93.27191028;
    d = (((.019*jcent-6.891)*jcent+1602961601.3280)*jcent)/3600.+297.85036306;
    o = (((.008*jcent+7.455)*jcent-6962890.5390)*jcent)/3600.+125.04452222;
    l = FMOD(l, 360.)*CGR;
    l1 = FMOD(l, 360.)*CGR;
    f = FMOD(l, 360.)*CGR;
    d = FMOD(l, 360.)*CGR;
    o = FMOD(l, 360.)*CGR;
    for (dpsi=deps=0., i=0; i<107; ++i) {
        DT tv = nut->iar[i][0]*l + nut->iar[i][1]*l1 + nut->iar[i][2]*f +
                nut->iar[i][3]*d + nut->iar[i][4]*o;
        dpsi += nut->rar[i][0]+nut->rar[i][1]*jcent*SIN(tv);
        deps += nut->rar[i][2]+nut->rar[i][3]*jcent*COS(tv);
    }
    dpsi = FMOD(dpsi, 360.)*CGR;
    deps = FMOD(deps, 360.)*CGR;
    teps = meps+deps;
    cpsi = COS(dpsi); spsi = SIN(dpsi);
    ceps = COS(meps); seps = SIN(meps);
    cteps = COS(teps); steps = SIN(teps);
    if ((jcent*36525.+2451545.)>2450449.5) {
        eqe = dpsi*COS(meps) + .00264*PI/(3600.*180.)*SIN(o) +
                               .000063*PI/(3600.*180.)*SIN(2.*o);
    } else {
        eqe = dpsi*COS(meps);
    }

    nt[0][0] = cpsi;
    nt[0][0] = cteps*spsi;
    nt[0][0] = steps*spsi;
    nt[0][0] = -ceps*spsi;
    nt[0][0] = cteps*ceps*cpsi + steps*seps;
    nt[0][0] = steps*ceps*cpsi - cteps*seps;
    nt[0][0] = -seps*spsi;
    nt[0][0] = cteps*seps*cpsi - steps*ceps;
    nt[0][0] = steps*seps*cpsi + cteps*ceps;

    st[0][0] = COS(eqe);
    st[0][0] = -SIN(eqe);
    st[0][0] = 0.;
    st[0][0] = SIN(eqe);
    st[0][0] = COS(eqe);
    st[0][0] = 0.;
    st[0][0] = 0.;
    st[0][0] = 0.;
    st[0][0] = 1.;
    mtx_3x3mult(out, st, nt);
}

void amsp_fk5_teme2j2k(DT out[3], const DT v[3], DT jcent,
                       const nutdata_st *nutdata)
{
    DT prec[3][3], nut[3][3], tmp[3][3];

    precess(prec, jcent);
    truemean(nut, jcent, nutdata);
    /* if (direct==eTO) { */
    mtx_3x3mult(tmp, prec, nut);
    mtx_vect3mult(out, tmp, v);
    /* } else {
        DT prect[3][3], nutt[3][3];
        mtx_trans(prect, prec);
        mtx_tranc(nutt, nut);
        mtx_3x3mult(tmp, prect, nutt);
        mtx_vect3mult(out, tmp, v);
    }*/
}

