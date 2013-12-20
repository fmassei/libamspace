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
#include "sgp4.h"
#include "precision.h"
#include "sgputils.h"
#include "date.h"

int sgp4_init(sat_st *sat, const gravconst_st *g)
{
    DT ss, qms4;            /* semi-constants, depend on perigee */
    DT cosio, sinio, theta2, x3thm1, beta,
       pinvsq, tsi, eta2, eeta, psi2, coef, coef1,
       x1mth2, theta4, x1m5th, xhdot1;
    DT tmp1, tmp2, tmp3;
    sat->c = sat->o;
    /* setup very common variables */
    cosio = COS(sat->o.i);
    sinio = SIN(sat->o.i);
    theta2 = P2(cosio);
    x3thm1 = 3.*theta2-1.;
    beta = SQRT(1-P2(sat->o.e));
    /* un-kozai the mean motion (get omm, sma and perigee) */
    utils_unkozai(sat, g);
    /* find out if sat is very near of not (as we'd need to use the linear
     * variation dropping some terms) */
    sat->nearth.isVeryNear = 0;
    if ((sat->sma*(1.-sat->o.e)/g->aE)<(220./g->rE+g->aE))
        sat->nearth.isVeryNear = 1;
    /* get ss and qms4 */
    utils_getDensityConstants(sat, g, &ss, &qms4);
    /* initialize stuff, mainly drag coefficients and derivatives */
    pinvsq = 1./(P2(sat->sma)*P4(beta));
    tsi = 1./(sat->sma-ss);
    sat->nearth.eta = sat->sma*sat->o.e*tsi;
    eta2 = P2(sat->nearth.eta);
    eeta = sat->o.e*sat->nearth.eta;
    psi2 = FABS(1.-eta2);
    coef = qms4*P4(tsi);
    coef1 = coef/POW(psi2, 3.5);
    sat->C[1] = coef1*sat->omm*(sat->sma*(1.+1.5*eta2+eeta*(4.+eta2))+.75*
        g->k2*tsi/psi2*x3thm1*(8.+3.*eta2*(8.+eta2)));
    sat->C[0]=sat->bstar*sat->C[1];
    sat->C[2] = 0.;
    if (sat->o.e>1.e-4)
        sat->C[2] = coef*tsi*g->a30cof*sat->omm*g->aE*sinio/sat->o.e;
    x1mth2 = 1.-theta2;
    sat->C[3] = 2.*sat->omm*coef1*sat->sma*P2(beta)*(sat->nearth.eta*
        (2.+.5*eta2)+sat->o.e*(.5+2.*eta2)-2.*g->k2*tsi/
        (sat->sma*psi2)*(-3.*x3thm1*(1.-2.*eeta+eta2*
        (1.5-.5*eeta))+.75*x1mth2*(2.*eta2-eeta*
        (1.+eta2))*COS(2.*sat->o.argp)));
    sat->C[4]=2.*coef1*sat->sma*P2(beta)*(1.+2.75*(eta2+eeta)+eeta*eta2);
    theta4 = theta2*theta2;
    tmp1 = 3.*g->k2*pinvsq*sat->omm;
    tmp2 = tmp1*g->k2*pinvsq;
    tmp3 = 1.25*g->k4*pinvsq*pinvsq*sat->omm;
    sat->c_d.ma = sat->omm+.5*tmp1*beta*x3thm1+.0625*tmp2*beta*
        (13.-78.*theta2+137.*theta4);
    x1m5th = 1.-5.*theta2;
    sat->c_d.argp = -.5*tmp1*x1m5th+.0625*tmp2*(7.-114.*theta2+
        395.*theta4)+tmp3*(3.-36.*theta2+49.*theta4);
    xhdot1 = -tmp1*cosio;
    sat->c_d.node = xhdot1+(.5*tmp2*(4.-19.*theta2)+2.*tmp3*(3.-
        7.*theta2))*cosio;
    sat->ncof = 3.5*P2(beta)*xhdot1*sat->C[0];
    sat->Tcof[0] = 1.5*sat->C[0];
    if (FABS(1.+cosio)>1.5e-12)
        sat->xcof = .125*g->a30cof*sinio*(3.+5.*cosio)/(1.+cosio);
    else
        sat->xcof = .125*g->a30cof*sinio*(3.+5.*cosio)/1.5e-12;
    sat->ycof = .25*g->a30cof*sinio;
    if (!sat->isDeep) {
        /* near earth, calculate drag */
        sat->nearth.dw = sat->bstar*sat->C[2]*COS(sat->o.argp);
        sat->nearth.dm = 0.;
        if (sat->o.e>1.e-4)
            sat->nearth.dm = -(2./3.)*coef*sat->bstar*g->aE/eeta;
        sat->nearth.dmo = POW((1.+sat->nearth.eta*COS(sat->o.ma)),3);
        if (!sat->nearth.isVeryNear) {
            /* non-linear version */
            DT tmp, c02;
            c02 = P2(sat->C[0]);
            sat->nearth.D[0] = 4.*sat->sma*tsi*c02;
            tmp = sat->nearth.D[0]*tsi*sat->C[0]/3.;
            sat->nearth.D[1] = (17.*sat->sma+ss)*tmp;
            sat->nearth.D[2] = .5*tmp*sat->sma*tsi*(221.*sat->sma+31.*ss)*
                               sat->C[0];
            sat->Tcof[1] = sat->nearth.D[0]+2.*c02;
            sat->Tcof[2] = .25*(3.*sat->nearth.D[1]+sat->C[0]*
                               (12.*sat->nearth.D[0]+10.*c02));
            sat->Tcof[3] = .2*(3.*sat->nearth.D[2]+12.*sat->C[0]*
                        sat->nearth.D[1]+6.*P2(sat->nearth.D[0])+15.*c02*(
                        2.*sat->nearth.D[0]+c02));
        }
    } else {
        /* deep space, call dpinit */
        dpinit(sat, g);
    }
    return 0;
}
int sgp4(sat_st *sat, const gravconst_st *g)
{
    DT cosio, sinio, beta,
       x3thm1, x1mth2, x7thm1;              /* common vars */
    DT dt, dt2;                             /* time tick */
    DT A, XL, AXN, XLL, AYNL, XLT, AYN,
       ecose, esine, PL,
       R, deltaR, deltaRF, betal, u;
    DT sinepw, cosepw;                      /* result of Kepler's equation */
    DT sinu, cosu, sin2u, cos2u;
    DT tmp, tmp1, tmp2, tmp3,
       tmpa, tmpe, tmpl;                    /* temp stuff */
    dt = sat->t;
    dt2 = P2(dt);
    /* very common stuff */
    cosio = COS(sat->o.i);
    sinio = SIN(sat->o.i);
    x1mth2 = 1.-P2(cosio);
    x3thm1 = 3.*P2(cosio)-1.;
    x7thm1 = 7.*P2(cosio)-1.;
    /* update for secular gravity and atmospheric drag */
    sat->c.ma = sat->o.ma+sat->c_d.ma*dt;
    sat->c.argp = sat->o.argp+sat->c_d.argp*dt;
    sat->c.node = (sat->o.node+sat->c_d.node*dt)+sat->ncof*dt2;
    tmpa = 1.-sat->C[0]*dt;
    tmpe = sat->bstar*sat->C[3]*dt;
    tmpl = sat->Tcof[0]*dt2;
    if (!sat->isDeep) {
        DT dt3, dt4, sinma;
        if (!sat->nearth.isVeryNear) {
            DT deltaArgp, deltaMM;
            deltaArgp = sat->nearth.dw*dt;
            deltaMM = sat->nearth.dm*(
                    POW((1.+sat->nearth.eta*COS(sat->c.ma)),3)-
                    sat->nearth.dmo);
            sat->c.ma += deltaArgp+deltaMM;
            sat->c.argp -= deltaArgp+deltaMM;
            dt3 = P3(dt);
            dt4 = P4(dt);
            tmpa -= sat->nearth.D[0]*dt2-
                      sat->nearth.D[1]*dt3-sat->nearth.D[2]*dt4;
            sinma = SIN(sat->o.ma);
            tmpe += sat->bstar*sat->C[4]*(SIN(sat->c.ma)-sinma);
            tmpl += sat->Tcof[1]*dt3+dt4*(sat->Tcof[2]+dt*sat->Tcof[3]);
        }
        A = sat->sma*P2(tmpa);
        sat->c.e = sat->o.e-tmpe;
        if (sat->c.e>=1.0 || sat->c.e<-0.001) {
            sat->error = SATERR_E;
            return -1;
        }
        XL = sat->c.ma+sat->c.argp+sat->c.node+sat->omm*tmpl;
        sat->c.i = sat->o.i;
    } else {
        sat->c.mm = sat->omm;
        dpsec(sat, &sat->c.ma, &sat->c.argp, &sat->c.node,
                &sat->c.e, &sat->c.i, &sat->c.mm);
        if (sat->c.mm<=0.) {
            sat->error = SATERR_MM;
            return -1;
        }
        A = POW(g->ke/sat->c.mm,2./3.)*P2(tmpa);
        sat->c.e -= tmpe;
        if (sat->c.e>=1.0 || sat->c.e<-0.001) {
            sat->error = SATERR_E;
            return -1;
        }
        sat->c.ma += sat->omm*tmpl;
        dpper(sat, g, &sat->c.e, &sat->c.i, &sat->c.argp, &sat->c.node,
                &sat->c.ma);
        if (sat->c.i<0.) {
            sat->c.i = -sat->c.i;
            sat->c.node += PI;
            sat->c.argp -= PI;
        }
        XL = sat->c.ma+sat->c.argp+sat->c.node;
        if (sat->c.e<.0 || sat->c.e>1.0) {
            sat->error = SATERR_E;
            return -1;
        }
    }
    if (sat->isDeep) {
        sinio = sin(sat->c.i);
        cosio = cos(sat->c.i);
        sat->ycof = .25*g->a30cof*sinio;
        if (FABS(1.+cosio)>1.5e-12)
            sat->xcof = .125*g->a30cof*sinio*(3.+5.*cosio)/(1.+cosio);
        else
            sat->xcof = .125*g->a30cof*sinio*(3.+5.*cosio)/1.5e-12;
        x1mth2 = 1.-P2(cosio);
        x3thm1 = 3.*P2(cosio)-1.;
        x7thm1 = 7.*P2(cosio)-1.;
    }
    beta = SQRT(1.-P2(sat->c.e));
    sat->c.mm = g->ke/POW(A, 1.5);
    /* long period periodics */
    AXN = sat->c.e*COS(sat->c.argp);
    tmp = 1./(A*P2(beta));
    XLL = tmp*sat->xcof*AXN;
    AYNL = tmp*sat->ycof;
    XLT = XL+XLL;
    AYN = sat->c.e*SIN(sat->c.argp)+AYNL;
    /* solve Kepler's equation */
    {
        DT u, epw, delta, it;
        u = FMOD(XLT-sat->c.node, 2.*PI);
        epw = u;
        delta = 1e100;
        it = 0;
        while ((FABS(delta)>=ZERO) && (it<10)) {
            sinepw = SIN(epw);
            cosepw = COS(epw);
            delta = (u-AYN*cosepw+AXN*sinepw-epw)/
                    (1.-cosepw*AXN-sinepw*AYN);
            if (FABS(delta)>=.95)
                delta = delta>0. ? .95 : -.95;
            epw += delta;
            ++it;
        }
    }
    ecose = AXN*cosepw+AYN*sinepw;
    esine = AXN*sinepw-AYN*cosepw;
    tmp = 1.-(P2(AXN)+P2(AYN));
    PL = A*tmp;
    if (PL<0.) {
        sat->error = SATERR_PL;
        return -1;
    }
    R = A*(1.-ecose);
    tmp1 = 1./R;
    deltaR = g->ke*SQRT(A)*esine*tmp1;
    deltaRF = g->ke*SQRT(PL)*tmp1;
    tmp2 = A*tmp1;
    betal = SQRT(tmp);
    tmp3 = 1./(1.+betal);
    cosu = tmp2*(cosepw-AXN+AYN*esine*tmp3);
    sinu = tmp2*(sinepw-AYN-AXN*esine*tmp3);
    u = ATAN2(sinu,cosu);
    sin2u = 2.*sinu*cosu;
    cos2u = 2.*cosu*cosu-1.;
    tmp = 1./PL;
    tmp1 = g->k2*tmp;
    tmp2 = tmp1*tmp;
    /* short periodics updates */
    R = R*(1.-1.5*tmp2*betal*x3thm1)+.5*tmp1*x1mth2*cos2u;
    u -= .25*tmp2*x7thm1*sin2u;
    sat->c.node += 1.5*tmp2*cosio*sin2u;
    sat->c.i += 1.5*tmp2*cosio*sinio*cos2u;
    deltaR -= sat->c.mm*tmp1*x1mth2*sin2u;
    deltaRF += sat->c.mm*tmp1*(x1mth2*cos2u+1.5*x3thm1);
    {
        DT xmx, xmy;
        DT ux, uy, uz, vx, vy, vz;              /* unit vector */
        DT sini, cosi, sinn, cosn;              /* just to read better (NOTE:
                                                 * sinu/cosu are overrided! */
        sinu = SIN(u);
        cosu = COS(u);
        sini = SIN(sat->c.i);
        cosi = COS(sat->c.i);
        sinn = SIN(sat->c.node);
        cosn = COS(sat->c.node);
        xmx = -sinn*cosi;
        xmy = cosn*cosi;
        /* get unit vectors */
        ux = xmx*sinu+cosn*cosu;
        uy = xmy*sinu+sinn*cosu;
        uz = sini*sinu;
        vx = xmx*cosu-cosn*sinu;
        vy = xmy*cosu-sinn*sinu;
        vz = sini*cosu;
        /* set position and velocity */
        sat->pos[0] = R*ux;
        sat->pos[1] = R*uy;
        sat->pos[2] = R*uz;
        sat->vel[0] = deltaR*ux+deltaRF*vx;
        sat->vel[1] = deltaR*uy+deltaRF*vy;
        sat->vel[2] = deltaR*uz+deltaRF*vz;
        VECT3_SCALAR(sat->pos, g->rE, sat->pos);
        VECT3_SCALAR(sat->vel, g->rE/60., sat->vel);
        if (R<1.0) {
            sat->error = SATERR_DECAY;
            return -1;
        }
    }
    return 0;
}

