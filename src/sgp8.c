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
#include "sgp8.h"
#include "precision.h"
#include "sgputils.h"
#include "date.h"

int sgp8_init(sat_st *sat, const gravconst_st *g)
{
    DT cosi, theta2, tthmun, e2, beta2, beta, B,
       po, pom2, sini, sing, cosg, tmp, theta4, unm5th, unmth2,
       pardt1, pardt2, pardt4,
       tsi, eta, eta2, psim2, alpha2, eeta, cos2g,
       D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16,
       D1DT, D2DT, D3DT, D4DT, D25, D17,
       B1, B2, B3, C0, C1, C4, C5, C4DT, C5DT,
       xndtn, D20, aldtal, tsdtts, etdt, psdtps,
       sin2g, C0DTC0, C1DTC1, xnddt, eddot, tsddts, etddt,
       C8, C9, D5DT, D18, D19, D23, D1DDT, xntrdt, tmnddt;
    const DT rho = .15696615;

    sat->c = sat->o;

    cosi = COS(sat->o.i);
    sini = SIN(sat->o.i);
    theta2 = P2(cosi);
    tthmun = 3.*theta2-1.;
    e2 = P2(sat->o.e);
    beta2 = 1.-e2;
    beta = SQRT(beta2);
    /* un-kozai the mean motion */
    utils_unkozai(sat, g);
    /* ballistic coefficient from B* drag term */
    B = 2.*sat->bstar/rho;
    /*INITIALIZATION*/
    sat->nearth.isVeryNear = 0;
    po = sat->sma*beta2;
    pom2 = 1./P2(po);
    sing = SIN(sat->o.argp);
    cosg = COS(sat->o.argp);
    theta4 = P2(theta2);
    unm5th = 1.-5.*theta2;
    unmth2 = 1.-theta2;
    pardt1 = 3.*g->k2*pom2*sat->omm;
    pardt2 = pardt1*g->k2*pom2;
    pardt4 = 1.25*g->k4*P2(pom2)*sat->omm;
    sat->mdt = .5*pardt1*beta*tthmun;
    sat->adt = -.5*pardt1*unm5th;
    sat->ndt = -pardt1*cosi;
    sat->c_d.ma = sat->omm+sat->mdt+
        .0625*pardt2*beta*(13.-78.*theta2+137.*theta4);
    sat->c_d.argp = sat->adt+
        .0625*pardt2*(7.-114.*theta2+395.*theta4)+pardt4*(3.-36.*
        theta2+49.*theta4);
    sat->c_d.node = sat->ndt+
        (.5*pardt2*(4.-19.*theta2)+2.*pardt4*(3.-7.*theta2))*cosi;
    tsi = 1./(po-g->s);
    eta = sat->o.e*g->s*tsi;
    eta2 = P2(eta);
    psim2 = FABS(1./(1.-eta2));
    alpha2 = 1.+e2;
    eeta = sat->o.e*eta;
    cos2g = 2.*P2(cosg)-1.;
    D5 = tsi*psim2;
    D1 = D5/po;
    D2 = 12.+eta2*(36.+4.5*eta2);
    D3 = eta2*(15.+2.5*eta2);
    D4 = eta*(5.+3.75*eta2);
    B1 = g->k2*tthmun;
    B2 = -g->k2*unmth2;
    B3 = g->a30cof*sini;
    C0 = .5*B*rho*g->qms4*sat->omm*sat->sma*P4(tsi)*POW(psim2,3.5)/SQRT(alpha2);
    C1 = 1.5*sat->omm*P2(alpha2)*C0;
    C4 = D1*D3*B2;
    C5 = D5*D4*B3;
    sat->c_d.mm = C1*(
        (2.+eta2*(3.+34.*e2)+5.*eeta*(4.+eta2)+8.5*e2)+
        D1*D2*B1+C4*cos2g+C5*sing);
    xndtn = sat->c_d.mm/sat->omm;
    /* if drag is very small the equations are truncated to linear variation
     * in mean motion and quadratic variation in mean anomaly */
    if (!sat->isDeep && FABS(xndtn*1440.)>=2.16E-3) {
        D6 = eta*(30.+22.5*eta2);
        D7 = eta*(5.+12.5*eta2);
        D8 = 1.+eta2*(6.75+eta2);
        C8 = D1*D7*B2;
        C9 = D5*D8*B3;
        sat->c_d.e = -C0*(
            eta*(4.+eta2+e2*(15.5+7.*eta2))+sat->o.e*(5.+15.*eta2)+
            D1*D6*B1 +
            C8*cos2g+C9*sing);
        D20 = .5*(2./3.)*xndtn;
        aldtal = sat->o.e*sat->c_d.e/alpha2;
        tsdtts = 2.*sat->sma*tsi*(D20*beta2+sat->o.e*sat->c_d.e);
        etdt = (sat->c_d.e+sat->o.e*tsdtts)*tsi*g->s;
        psdtps = -eta*etdt*psim2;
        sin2g = 2.*sing*cosg;
        C0DTC0 = D20+4.*tsdtts-aldtal-7.*psdtps;
        C1DTC1 = xndtn+4.*aldtal+C0DTC0;
        D9 = eta*(6.+68.*e2)+sat->o.e*(20.+15.*eta2);
        D10 = 5.*eta*(4.+eta2)+sat->o.e*(17.+68.*eta2);
        D11 = eta*(72.+18.*eta2);
        D12 = eta*(30.+10.*eta2);
        D13 = 5.+11.25*eta2;
        D14 = tsdtts-2.*psdtps;
        D15 = 2.*(D20+sat->o.e*sat->c_d.e/beta2);
        D1DT = D1*(D14+D15);
        D2DT = etdt*D11;
        D3DT = etdt*D12;
        D4DT = etdt*D13;
        D5DT = D5*D14;
        C4DT = B2*(D1DT*D3+D1*D3DT);
        C5DT = B3*(D5DT*D4+D5*D4DT);
        D16 = D9*etdt+D10*sat->c_d.e +
            B1*(D1DT*D2+D1*D2DT) +
            C4DT*cos2g+C5DT*sing+sat->adt*(C5*cosg-2.*C4*sin2g);
        xnddt = C1DTC1*sat->c_d.mm+C1*D16;
        eddot = C0DTC0*sat->c_d.e-C0*(
            (4.+3.*eta2+30.*eeta+e2*(15.5+21.*eta2))*etdt+(5.+15.*eta2
            +eeta*(31.+14.*eta2))*sat->c_d.e +
            B1*(D1DT*D6+D1*etdt*(30.+67.5*eta2)) +
            B2*(D1DT*D7+D1*etdt*(5.+37.5*eta2))*cos2g+
            B3*(D5DT*D8+D5*etdt*eta*(13.5+4.*eta2))*sing+sat->adt*(C9*
            cosg-2.*C8*sin2g));
        D25 = P2(sat->c_d.e);
        D17 = xnddt/sat->omm-P2(xndtn);
        tsddts = 2.*tsdtts*(tsdtts-D20)+sat->sma*tsi*((2./3.)*beta2*D17-4.*D20*
            sat->o.e*sat->c_d.e+2.*(D25+sat->o.e*eddot));
        etddt =(eddot+2.*sat->c_d.e*tsdtts)*tsi*g->s+tsddts*eta;
        D18 = tsddts-P2(tsdtts);
        D19 = -P2(psdtps)/eta2-eta*etddt*psim2-P2(psdtps);
        D23 = P2(etdt);
        D1DDT = D1DT*(D14+D15)+D1*(D18-2.*D19+(2./3.)*D17+2.*(alpha2*D25
            /beta2+sat->o.e*eddot)/beta2);
        xntrdt = sat->c_d.mm*(2.*(2./3.)*D17+3.*
            (D25+sat->o.e*eddot)/alpha2-6.*P2(aldtal) +
            4.*D18-7.*D19 )+
            C1DTC1*xnddt+C1*(C1DTC1*D16+
            D9*etddt+D10*eddot+D23*(6.+30.*eeta+68.*e2)+
            etdt*sat->c_d.e*(40.+30.*
            eta2+272.*eeta)+D25*(17.+68.*eta2) +
            B1*(D1DDT*D2+2.*D1DT*D2DT+D1*(etddt*D11+D23*(72.+54.*eta2))) +
            B2*(D1DDT*D3+2.*D1DT*D3DT+D1*(etddt*D12+D23*(30.+30.*eta2))) *
            cos2g+
            B3*((D5DT*D14+D5*(D18-2.*D19)) *
            D4+2.*D4DT*D5DT+D5*(etddt*D13+22.5*eta*D23)) *sing+sat->adt*
            ((7.*D20+4.*sat->o.e*sat->c_d.e/beta2)*
            (C5*cosg-2.*C4*sin2g)
            +((2.*C5DT*cosg-4.*C4DT*sin2g)-sat->adt*(C5*sing+4.*
            C4*cos2g))));
        tmnddt = xnddt*1.E9;
        tmp = P2(tmnddt)-sat->c_d.mm*1.E18*xntrdt;
        sat->pp = (tmp+P2(tmnddt))/tmp;
        sat->gamma = -xntrdt/(xnddt*(sat->pp-2.));
        sat->d_omm = sat->c_d.mm/(sat->pp*sat->gamma);
        sat->qq = 1.-eddot/(sat->c_d.e*sat->gamma);
        sat->ed = sat->c_d.e/(sat->qq*sat->gamma);
        sat->ovgpp = 1./(sat->gamma*(sat->pp+1.));
    } else {
        sat->nearth.isVeryNear = 1;
        sat->c_d.e = -(2./3.)*xndtn*(1.-sat->o.e);
    }
    if (sat->isDeep) {
        dpinit(sat, g);
    }
    return 0;
}

int sgp8(sat_st *sat, const gravconst_st *g)
{
    DT cosi, theta2, tthmun, sini, sinio2, cosio2, unm5th, unmth2;
    DT Z1, Z7,
       sine, cose, ZC5, am, beta2m, sinargp, cosargp, axnm, aynm,
       pm, G1, G2, G3, beta, G4, G5, SNF, CSF, FM, sinu, cosu, SN2F2G, CS2F2G,
       ecosf, G10, RM, AOVR, G13, G14, DR, DIWC, DI, SNI2DU, lambda, Y4, Y5,
       sini2;
    DT tmp;
    DT dt;

    dt = sat->t;

    cosi = COS(sat->o.i);
    sini = SIN(sat->o.i);
    sinio2 = SIN(.5*sat->o.i);
    cosio2 = COS(.5*sat->o.i);
    theta2 = P2(cosi);
    tthmun = 3.*theta2-1.;
    unm5th = 1.-5.*theta2;
    unmth2 = 1.-theta2;
    /* update for secular gravity and atmospheric drag */
    if (!sat->isDeep) {
        sat->c.ma = sat->o.ma+sat->c_d.ma*dt;
        sat->c.argp = sat->o.argp+sat->c_d.argp*dt;
        sat->c.node = sat->o.node+sat->c_d.node*dt;
        sat->c.ma = FMOD(sat->c.ma, 2.*PI);
        if (!sat->nearth.isVeryNear) {
            DT tmp1;
            tmp = 1.-sat->gamma*dt;
            tmp1 = POW(tmp, sat->pp);
            sat->c.mm = sat->omm+sat->d_omm*(1.-tmp1);
            sat->c.e = sat->o.e+sat->ed*(1.-POW(tmp,sat->qq));
            Z1 = sat->d_omm*(dt+sat->ovgpp*(tmp*tmp1-1.));
        } else {
            sat->c.mm = sat->omm+sat->c_d.mm*dt;
            sat->c.e = sat->o.e+sat->c_d.e*dt;
            Z1 = .5*sat->c_d.mm*P2(dt);
        }
        Z7 = 3.5*(2./3.)*Z1/sat->omm;
        sat->c.ma += Z1+Z7*sat->mdt;
        sat->c.argp += Z7*sat->adt;
        sat->c.node += Z7*sat->ndt;
        sat->c.ma = FMOD(sat->c.ma, 2.*PI);
        if (sat->c.e>=1. || sat->c.e<-0.001) {
            sat->error = SATERR_E;
            return -1;
        }
    } else {
        Z1 = .5*sat->c_d.mm*P2(dt);
        Z7 = 3.5*(2./3.)*Z1/sat->omm;
        sat->c.ma = sat->o.ma+sat->c_d.ma*dt;
        sat->c.argp = sat->o.argp+sat->c_d.argp*dt+Z7*sat->adt;
        sat->c.node = sat->o.node+sat->c_d.node*dt+Z7*sat->ndt;
        sat->c.mm = sat->omm;
        dpsec(sat, &sat->c.ma, &sat->c.argp, &sat->c.node, &sat->c.e,
                   &sat->c.i, &sat->c.mm);
        if (sat->c.mm<=0.) {
            sat->error = SATERR_MM;
            return -1;
        }
        sat->c.mm += sat->c_d.mm*dt;
        sat->c.e += sat->c_d.e*dt;
        sat->c.ma += Z1+Z7*sat->mdt;
        if (sat->c.e>=1. || sat->c.e<-0.001) {
            sat->error = SATERR_E;
            return -1;
        }
        dpper(sat, g, &sat->c.e, &sat->c.i, &sat->c.argp, &sat->c.node,
                      &sat->c.ma);
        sat->c.ma = FMOD(sat->c.ma, 2.*PI);
        if (sat->c.i<0.) {
            sat->c.i = -sat->c.i;
            sat->c.node += PI;
            sat->c.argp -= PI;
        }
        if (sat->c.e>=1. || sat->c.e<-0.001) {
            sat->error = SATERR_E;
            return -1;
        }
    }
    /* solve Kepler's equation */
    {
        DT delta = 1e100, u, epw;
        int it = 0;
        u = sat->c.ma+sat->c.e*SIN(sat->c.ma)*(1+sat->c.e*COS(sat->c.ma));
        epw = u;
        while ((FABS(delta)>=ZERO)&&(it<10)) {
            sine = SIN(epw);
            cose = COS(epw);
            ZC5 = 1./(1-sat->c.e*cose);
            delta = (sat->c.ma+sat->c.e*sine-epw)*ZC5;
            if (FABS(delta)>=.95)
                delta = delta>0. ? .95 : -95;
            epw += delta;
            ++it;
        }
    }
    /* short period preliminary quantities */
    am = POW(g->ke/sat->c.mm,(2./3.));
    beta2m = 1.-P2(sat->c.e);
    sinargp = SIN(sat->c.argp);
    cosargp = COS(sat->c.argp);
    axnm = sat->c.e*cosargp;
    aynm = sat->c.e*sinargp;
    pm = am*beta2m;
    G1 = 1./pm;
    G2 = .5*g->k2*G1;
    G3 = G2*G1;
    beta = SQRT(beta2m);
    G4 = .25*g->a30cof*sini;
    G5 = .25*g->a30cof*G1;
    SNF = beta*sine*ZC5;
    CSF = (cose-sat->c.e)*ZC5;
    FM = ATAN2(SNF,CSF);
    sinu = SNF*cosargp+CSF*sinargp;
    cosu = CSF*cosargp-SNF*sinargp;
    SN2F2G = 2.*sinu*cosu;
    CS2F2G = 2.*P2(cosu)-1.;
    ecosf = sat->c.e*CSF;
    G10 = FM-sat->c.ma+sat->c.e*SNF;
    RM = pm/(1.+ecosf);
    AOVR = am/RM;
    G13 = sat->c.mm*AOVR;
    G14 = -G13*AOVR;
    DR = G2*(unmth2*CS2F2G-3.*tthmun)-G4*sinu;
    DIWC = 3.*G3*sini*CS2F2G-G5*aynm;
    DI = DIWC*cosi;
    sini2 = SIN(.5*sat->c.i);
    /* update for short period periodics */
    SNI2DU = sinio2*(
        G3*(.5*(1.-7.*theta2)*SN2F2G-3.*unm5th*G10)-G5*sini*cosu*(2.+
        ecosf))-.5*G5*theta2*axnm/cosio2;
    lambda = sat->c.argp+sat->c.node+FM+G3*(.5*(1.+6.*cosi-7.*theta2)*SN2F2G-3.*
        (unm5th+2.*cosi)*G10)+G5*sini*(cosi*axnm/(1.+cosi)-(2.
        +ecosf)*cosu);
    if (!sat->isDeep) {
        Y4 = sinio2*sinu+cosu*SNI2DU+.5*sinu*cosio2*DI;
        Y5 = sinio2*cosu-sinu*SNI2DU+.5*cosu*cosio2*DI;
    } else {
        Y4 = sini2*sinu+cosu*SNI2DU+.5*sinu*cosio2*DI;
        Y5 = sini2*cosu-sinu*SNI2DU+.5*cosu*cosio2*DI;
    }
    {
        DT ux, uy, uz, vx, vy, vz;
        DT snlamb, cslamb, rdot, rvdot, r;

        r = RM+DR;
        rdot = sat->c.mm*am*sat->c.e*SNF/beta+G14*(2.*G2*unmth2*SN2F2G+G4*cosu);
        rvdot = sat->c.mm*P2(am)*beta/RM+G14*DR+am*G13*sini*DIWC;
        /* get unit vectors */
        snlamb = SIN(lambda);
        cslamb = COS(lambda);
        tmp = 2.*(Y5*snlamb-Y4*cslamb);
        ux = Y4*tmp+cslamb;
        vx = Y5*tmp-snlamb;
        tmp = 2.*(Y5*cslamb+Y4*snlamb);
        uy = -Y4*tmp+snlamb;
        vy = -Y5*tmp+cslamb;
        tmp = 2.*SQRT(1.-P2(Y4)-P2(Y5));
        uz = Y4*tmp;
        vz = Y5*tmp;
        /* set position and velocity */
        sat->pos[0] = r*ux;
        sat->pos[1] = r*uy;
        sat->pos[2] = r*uz;
        sat->vel[0] = rdot*ux+rvdot*vx;
        sat->vel[1] = rdot*uy+rvdot*vy;
        sat->vel[2] = rdot*uz+rvdot*vz;
        VECT3_SCALAR(sat->pos, g->rE, sat->pos);
        VECT3_SCALAR(sat->vel, g->rE/60., sat->vel);
    }
    return 0;
}

