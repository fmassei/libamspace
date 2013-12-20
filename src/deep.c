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
#include "deep.h"
#include "sat.h"
#include "date.h"

#define ZCOSIS  .91744867
#define ZSINIS  .39785416
#define ZCOSGS  .1945905
#define ZSINGS  -.98088458

#define Q22     1.7891679E-6
#define Q31     2.1460748E-6
#define Q33     2.2123015E-7
#define G22     5.7686396
#define G32     0.95240898
#define G44     1.8014998
#define G52     1.0508330
#define G54     4.4108898
#define ROOT22  1.7891679E-6
#define ROOT32  3.7393792E-7
#define ROOT44  7.3636953E-9
#define ROOT52  1.1428639E-7
#define ROOT54  2.1765803E-9
#define THDT    4.3752691E-3

#define FASX2   .13130908
#define FASX4   2.8843198
#define FASX6   .37448087

#define STEP    720.

void dpinit(sat_st *sat, const gravconst_st *g)
{
    DT F[6][5][4];      /* (F_lmp) inclination function */
    DT G[6][4][4];      /* (G_lpq) eccentricity function */
    DT A1, A2, A3, A4, A5, A6, A7, A8, A9, A10,
       X1, X2, X3, X4, X5, X6, X7, X8,
       Z11, Z12, Z13, Z21, Z22, Z23, Z31, Z32, Z33, Z1, Z2, Z3,
       S1, S2, S3, S4, S5, S6, S7,
       ZX, ZY, ZN, ZE, SE, SI, SL, SGH, SH; /* all these vars have no physical
                                             * meaning.. */
    DT temp1, temp;
    DT gamma, bfact, sini, cosi, theta, beta;
    DT e, e2, e3;       /* just to read better: e, e^2 and e^3 */
    DT sinnode, cosnode, sinargp, cosargp, sini2, aqnv;
    DT zcosil, zsinil, zcoshl, zsinhl, zcosg, zsing, zcosh, zsinh,
       zcosgl, zsingl, zsini, zcosi;
    DT C, CC;
    DT xpidot, xnodce, stem, ctem;
    DT day, doingSun;
    sat->deep.gst = amsp_date_gstime(sat->epochJD);
    
    e = sat->o.e;
    e2 = P2(e);
    e3 = P3(e);
    aqnv = 1./sat->sma;
    xpidot = sat->c_d.argp+sat->c_d.node;
    sini = SIN(sat->o.i);
    cosi = COS(sat->o.i);
    theta = P2(cosi);
    sinnode = SIN(sat->o.node);
    cosnode = COS(sat->o.node);
    sinargp = SIN(sat->o.argp);
    cosargp = COS(sat->o.argp);
    beta = SQRT(1-P2(sat->o.e));
    /*INITIALIZE LUNAR SOLAR TERMS*/
    day = sat->epochJD-2433281.5+18261.5;
    xnodce = 4.5236020-9.2422029e-4*day;
    stem = SIN(xnodce);
    ctem = COS(xnodce);
    zcosil=.91375164-.03568096*ctem;
    zsinil=SQRT (1.-P2(zcosil));
    zsinhl= .089683511*stem/zsinil;
    zcoshl=SQRT (1.-P2(zsinhl));
    C=4.7199672+.22997150*day;
    gamma=5.8351514+.0019443680*day;
    sat->deep.moon.zmo = FMOD(C-gamma, 2.*PI);
    ZX= .39785416*stem/zsinil;
    ZY= zcoshl*ctem+0.91744867*zsinhl*stem;
    ZX=ATAN2(ZX,ZY);
    ZX=gamma+ZX-xnodce;
    zcosgl=COS (ZX);
    zsingl=SIN (ZX);
    sat->deep.sun.zmo=FMOD(6.2565837+.017201977*day, 2.*PI);
    /*DO SOLAR TERMS*/
    zcosg=ZCOSGS;
    zsing=ZSINGS;
    zcosi=ZCOSIS;
    zsini=ZSINIS;
    zcosh=cosnode;
    zsinh=sinnode;
    CC=g->cS;
    ZN=g->nS;
    ZE=g->eS;
    doingSun = 1;
    gravloop:
    A1=zcosg*zcosh+zsing*zcosi*zsinh;
    A3=-zsing*zcosh+zcosg*zcosi*zsinh;
    A7=-zcosg*zsinh+zsing*zcosi*zcosh;
    A8=zsing*zsini;
    A9=zsing*zsinh+zcosg*zcosi*zcosh;
    A10=zcosg*zsini;
    A2= cosi*A7+ sini*A8;
    A4= cosi*A9+ sini*A10;
    A5=- sini*A7+ cosi*A8;
    A6=- sini*A9+ cosi*A10;
    
    X1=A1*cosargp+A2*sinargp;
    X2=A3*cosargp+A4*sinargp;
    X3=-A1*sinargp+A2*cosargp;
    X4=-A3*sinargp+A4*cosargp;
    X5=A5*sinargp;
    X6=A6*sinargp;
    X7=A5*cosargp;
    X8=A6*cosargp;

    Z31=12.*X1*X1-3.*X3*X3;
    Z32=24.*X1*X2-6.*X3*X4;
    Z33=12.*X2*X2-3.*X4*X4;
    Z1=3.*(A1*A1+A2*A2)+Z31*e2;
    Z2=6.*(A1*A3+A2*A4)+Z32*e2;
    Z3=3.*(A3*A3+A4*A4)+Z33*e2;
    Z11=-6.*A1*A5+e2 *(-24.*X1*X7-6.*X3*X5);
    Z12=-6.*(A1*A6+A3*A5)+e2 *(-24.*(X2*X7+X1*X8)-6.*(X3*X6+X4*X5));
    Z13=-6.*A3*A6+e2 *(-24.*X2*X8-6.*X4*X6);
    Z21=6.*A2*A5+e2 *(24.*X1*X5-6.*X3*X7);
    Z22=6.*(A4*A5+A2*A6)+e2 *(24.*(X2*X5+X1*X6)-6.*(X4*X7+X3*X8));
    Z23=6.*A4*A6+e2 *(24.*X2*X6-6.*X4*X8);
    Z1=Z1+Z1+P2(beta)*Z31;
    Z2=Z2+Z2+P2(beta)*Z32;
    Z3=Z3+Z3+P2(beta)*Z33;
    S3=CC*(1./sat->omm);
    S2=-.5*S3/beta;
    S4=S3*beta;
    S1=-15.*e*S4;
    S5=X1*X3+X2*X4;
    S6=X2*X3+X1*X4;
    S7=X2*X4-X1*X3;
    SE=S1*ZN*S5;
    SI=S2*ZN*(Z11+Z13);
    SL=-ZN*S3*(Z1+Z3-14.-6.*e2);
    SGH=S4*ZN*(Z31+Z33-6.);
    SH=-ZN*S2*(Z21+Z23);
    if (sat->o.i<5.2359877E-2 || sat->o.i>(PI-5.2359877E-2)) SH=0.0;
    sat->deep.moon.e[0]=2.*S1*S6;
    sat->deep.moon.e[1]=2.*S1*S7;
    sat->deep.moon.i[0]=2.*S2*Z12;
    sat->deep.moon.i[1]=2.*S2*(Z13-Z11);
    sat->deep.moon.l[0]=-2.*S3*Z2;
    sat->deep.moon.l[1]=-2.*S3*(Z3-Z1);
    sat->deep.moon.l[2]=-2.*S3*(-21.-9.*e2)*ZE;
    sat->deep.moon.gh[0]=2.*S4*Z32;
    sat->deep.moon.gh[1]=2.*S4*(Z33-Z31);
    sat->deep.moon.gh[2]=-18.*S4*ZE;
    sat->deep.moon.h[0]=-2.*S2*Z22;
    sat->deep.moon.h[1]=-2.*S2*(Z23-Z21);
    if (doingSun) goto l30; else goto l40;
    /*DO LUNAR TERMS*/
    l30:
    sat->deep.de.e = SE;
    sat->deep.de.i = SI;
    sat->deep.de.ma = SL;
    sat->deep.de.node=SH/sini;
    sat->deep.de.argp=SGH-cosi*sat->deep.de.node;
    sat->deep.sun.e[0]=sat->deep.moon.e[0];
    sat->deep.sun.i[0]=sat->deep.moon.i[0];
    sat->deep.sun.l[0]=sat->deep.moon.l[0];
    sat->deep.sun.gh[0]=sat->deep.moon.gh[0];
    sat->deep.sun.h[0]=sat->deep.moon.h[0];
    sat->deep.sun.e[1]=sat->deep.moon.e[1];
    sat->deep.sun.i[1]=sat->deep.moon.i[1];
    sat->deep.sun.l[1]=sat->deep.moon.l[1];
    sat->deep.sun.gh[1]=sat->deep.moon.gh[1];
    sat->deep.sun.h[1]=sat->deep.moon.h[1];
    sat->deep.sun.l[2]=sat->deep.moon.l[2];
    sat->deep.sun.gh[2]=sat->deep.moon.gh[2];
    zcosg=zcosgl;
    zsing=zsingl;
    zcosi=zcosil;
    zsini=zsinil;
    zcosh=zcoshl*cosnode+zsinhl*sinnode;
    zsinh=sinnode*zcoshl-cosnode*zsinhl;
    ZN=g->nL;
    CC=g->cL;
    ZE=g->eL;
    doingSun = 0;
    goto gravloop;
    l40:
    sat->deep.de.e += SE;
    sat->deep.de.i += SI;
    sat->deep.de.ma += SL;
    sat->deep.de.argp += SGH-cosi/sini*SH;
    sat->deep.de.node += SH/sini;
    /*GEOPOTENTIAL RESONANCE INITIALIZATION FOR 12 HOUR ORBITS*/
    sat->deep.res = RESONANCE_NONE;
    if (sat->omm<.0052359877 && sat->omm>.0034906585) {
        sat->deep.res = RESONANCE_ONEDAY;
    }
    if (sat->omm>=8.26E-3 && sat->omm<=9.24E-3 && e>=.5) {
        sat->deep.res = RESONANCE_HALFDAY;
    }
    if (sat->deep.res==RESONANCE_NONE)
        return;
    if (sat->deep.res==RESONANCE_ONEDAY) {
        G[2][0][0] = 1.0+e2*(-2.5+.8125*e2);
        G[3][1][0] = 1.0+2.0*e2;
        G[3][0][0] = 1.0+e2*(-6.0+6.60937*e2);
        F[2][2][0] = .75*(1.+cosi)*(1.+cosi);
        F[3][1][1] = .9375*P2(sini)*(1.+3.*cosi)-.75*(1.+cosi);
        F[3][3][0] = 1.+cosi;
        F[3][3][0] = 1.875*P3(F[3][3][0]);
        sat->deep.el[0] = 3.*P2(sat->omm)*P2(aqnv);
        sat->deep.el[1] = 2.*sat->deep.el[0]*F[2][2][0]*G[2][0][0]*Q22;
        sat->deep.el[2] = 3.*sat->deep.el[0]*F[3][3][0]*G[3][0][0]*Q33*aqnv;
        sat->deep.el[0] = sat->deep.el[0]*F[3][1][1]*G[3][1][0]*Q31*aqnv;
        sat->deep.lam = FMOD(sat->o.ma + sat->o.node + sat->o.argp
                            -FMOD(sat->deep.gst, 2.*PI), 2.*PI);
        bfact = sat->c_d.ma + xpidot - THDT +
                sat->deep.de.ma + sat->deep.de.argp + sat->deep.de.node;
    } else { /* RESONANCE_HALFDAY */
        G[2][0][1] = -.306-(e-.64)*.440;
        if (e<=.65) {
            G[2][1][1] = 3.616-13.247*e+16.290*e2;
            G[3][1][0] = -19.302+117.390*e-228.419*e2+156.591*e3;
            G[3][2][2] = -18.9068+109.7927*e-214.6334*e2+146.5816*e3;
            G[4][1][0] = -41.122+242.694*e-471.094*e2+313.953*e3;
            G[4][2][2] = -146.407+841.880*e-1629.014*e2+1083.435*e3;
            G[5][2][0] = -532.114+3017.977*e-5740*e2+3708.276*e3;
        } else {
            G[2][1][1] = -72.099+331.819*e-508.738*e2+266.724*e3;
            G[3][1][0] = -346.844+1582.851*e-2415.925*e2+1246.113*e3;
            G[3][2][2] = -342.585+1554.908*e-2366.899*e2+1215.972*e3;
            G[4][1][0] = -1052.797+4758.686*e-7193.992*e2+3651.957*e3;
            G[4][2][2] = -3581.69+16178.11*e-24462.77*e2+12422.52*e3;
            if (e>.715)
                G[5][2][0] = -5149.66+29936.92*e-54087.36*e2+31324.56*e3;
            else
                G[5][2][0] = 1464.74-4664.75*e+3763.64*e2;
        }
        if (e<.7) {
            G[5][3][3] = -919.2277+4988.61*e-9064.77*e2+5542.21*e3;
            G[5][2][1] = -822.71072+4568.6173*e-8491.4146*e2+5337.524*e3;
            G[5][3][2] = -853.666+4690.25*e-8624.77*e2+5341.4*e3;
        } else {
            G[5][3][3] = -37995.78+161616.52*e-229838.2*e2+109377.94*e3;
            G[5][2][1] = -51752.104+218913.95*e-309468.16*e2+146349.42*e3;
            G[5][3][2] = -40023.88+170470.89*e-242699.48*e2+115605.82*e3;
        }
        sini2=P2(sini);
        F[2][2][0] = .75*(1.+2.*cosi+theta);
        F[2][2][1] = 1.5*sini2;
        F[3][2][1] = 1.875*sini*(1.-2.*cosi-3.*theta);
        F[3][2][2] = -1.875*sini*(1.+2.*cosi-3.*theta);
        F[4][4][1] = 35.*sini2*F[2][2][0];
        F[4][4][2] = 39.3750*sini2*sini2;
        F[5][2][2] = 9.84375*sini*(sini2*(1.-2.*cosi-5.*theta)
            +.33333333*(-2.+4.*cosi+6.*theta));
        F[5][2][3] = sini*(4.92187512*sini2*(-2.-4.*cosi+10.*theta)
            +6.56250012*(1.+2.*cosi-3.*theta));
        F[5][4][2] = 29.53125*sini*(2.-8.*cosi+theta*(-12.+8.*cosi
            +10.*theta));
        F[5][4][3] = 29.53125*sini*(-2.-8.*cosi+theta*(12.+8.*cosi-10.*theta));
        temp1 = 3.*P2(sat->omm)*P2(aqnv);
        temp = temp1*ROOT22;
        sat->deep.q[2][2][0][1] = temp*F[2][2][0]*G[2][0][1];
        sat->deep.q[2][2][1][1] = temp*F[2][2][1]*G[2][1][1];
        temp1 = temp1*aqnv;
        temp = temp1*ROOT32;
        sat->deep.q[3][2][1][0] = temp*F[3][2][1]*G[3][1][0];
        sat->deep.q[3][2][2][2] = temp*F[3][2][2]*G[3][2][2];
        temp1 = temp1*aqnv;
        temp = 2.*temp1*ROOT44;
        sat->deep.q[4][4][1][0] = temp*F[4][4][1]*G[4][1][0];
        sat->deep.q[4][4][2][2] = temp*F[4][4][2]*G[4][2][2];
        temp1 = temp1*aqnv;
        temp = temp1*ROOT52;
        sat->deep.q[5][2][2][0] = temp*F[5][2][2]*G[5][2][0];
        sat->deep.q[5][2][3][2] = temp*F[5][2][3]*G[5][3][2];
        temp = 2.*temp1*ROOT54;
        sat->deep.q[5][4][2][1] = temp*F[5][4][2]*G[5][2][1];
        sat->deep.q[5][4][3][3] = temp*F[5][4][3]*G[5][3][3];
        sat->deep.lam = FMOD(sat->o.ma+2.*sat->o.node-
                             2.*FMOD(sat->deep.gst, 2.*PI), 2.*PI);
        bfact = sat->c_d.ma + sat->deep.de.ma +
                2.*(sat->c_d.node+sat->deep.de.node-THDT);
    }
    sat->deep.fact = bfact-sat->omm;
}

void dpsec(sat_st *sat,
            DT *XLL, DT *OMGASM, DT *XNODES, DT *EM, DT *XINC, DT *XN)
{
    DT FT, XNDOT, XNDDT, XL, XLDOT;
    DT XLI, XNI, ATIME, STEP2;
    DT delta, temp;
    DT dt;
    int done=0;
    STEP2 = 259200;
    ATIME = 0;
    dt = sat->t;
    *XLL=*XLL+sat->deep.de.ma*dt;
    *OMGASM=*OMGASM+sat->deep.de.argp*dt;
    *XNODES=*XNODES+sat->deep.de.node*dt;
    *EM=sat->o.e+sat->deep.de.e*dt;
    *XINC=sat->o.i+sat->deep.de.i*dt;

    FT = 0.;
    if (sat->deep.res==RESONANCE_NONE) return;
    if ((ATIME==0.) || (dt*ATIME<=0.) || (FABS(dt)<FABS(ATIME))) {
        ATIME = 0.;
        XNI = sat->omm;
        XLI = sat->deep.lam;    
    }
    delta = (dt>=0.)?STEP:-STEP;
    while (!done) {
        if (sat->deep.res==RESONANCE_ONEDAY) {
            XNDOT=sat->deep.el[0]*SIN (XLI-FASX2)+sat->deep.el[1]*SIN (2.*(XLI-FASX4))
                +sat->deep.el[2]*SIN (3.*(XLI-FASX6));
            XNDDT = sat->deep.el[0]*COS(XLI-FASX2)
                +2.*sat->deep.el[1]*COS(2.*(XLI-FASX4))
                +3.*sat->deep.el[2]*COS(3.*(XLI-FASX6));
        } else {
            DT XOMI, X2OMI, X2LI;
            XOMI = sat->o.argp+sat->c_d.argp*ATIME;
            X2OMI = XOMI+XOMI;
            X2LI = XLI+XLI;
            XNDOT = sat->deep.q[2][2][0][1]*SIN(X2OMI+XLI-G22)
                +sat->deep.q[2][2][1][1]*SIN(XLI-G22)
                +sat->deep.q[3][2][1][0]*SIN(XOMI+XLI-G32)
                +sat->deep.q[3][2][2][2]*SIN(-XOMI+XLI-G32)
                +sat->deep.q[4][4][1][0]*SIN(X2OMI+X2LI-G44)
                +sat->deep.q[4][4][2][2]*SIN(X2LI-G44)
                +sat->deep.q[5][2][2][0]*SIN(XOMI+XLI-G52)
                +sat->deep.q[5][2][3][2]*SIN(-XOMI+XLI-G52)
                +sat->deep.q[5][4][2][1]*SIN(XOMI+X2LI-G54)
                +sat->deep.q[5][4][3][3]*SIN(-XOMI+X2LI-G54);
            XNDDT = sat->deep.q[2][2][0][1]*COS(X2OMI+XLI-G22)
                +sat->deep.q[2][2][1][1]*COS(XLI-G22)
                +sat->deep.q[3][2][1][0]*COS(XOMI+XLI-G32)
                +sat->deep.q[3][2][2][2]*COS(-XOMI+XLI-G32)
                +sat->deep.q[5][2][2][0]*COS(XOMI+XLI-G52)
                +sat->deep.q[5][2][3][2]*COS(-XOMI+XLI-G52)
                +2.*(sat->deep.q[4][4][1][0]*COS(X2OMI+X2LI-G44)
                +sat->deep.q[4][4][2][2]*COS(X2LI-G44)
                +sat->deep.q[5][4][2][1]*COS(XOMI+X2LI-G54)
                +sat->deep.q[5][4][3][3]*COS(-XOMI+X2LI-G54));
        }
        XLDOT=XNI+sat->deep.fact;
        XNDDT = XNDDT*XLDOT;
        if (FABS(dt-ATIME)<STEP) {
            FT = dt-ATIME;
            done=1;
        }
        if (!done) {
            XLI = XLI+XLDOT*delta+XNDOT*STEP2;
            XNI = XNI+XNDOT*delta+XNDDT*STEP2;
            ATIME=ATIME+delta;
        }
    }
    *XN = XNI+XNDOT*FT+XNDDT*FT*FT*0.5;
    XL = XLI+XLDOT*FT+XNDOT*FT*FT*0.5;
    temp = FMOD(sat->deep.gst+dt*THDT, 2.*PI);
    if (sat->deep.res==RESONANCE_HALFDAY) {
        *XLL = XL-2.*(*XNODES)+2.*temp;
    } else {
        *XLL = XL-*XNODES-*OMGASM+temp;
    }
}

void dpper(sat_st *sat, const gravconst_st *g,
            DT *EM, DT *XINC, DT *OMGASM, DT *XNODES, DT *XLL)
{
    DT sinis, cosis;
    DT ZM, ZF, f2, f3,
       SINZF,
       PE, PINC, PL, PGH, PH,
       SINOK, COSOK, ALFDP, BETDP, DALF, DBET, XLS, DLS;
    DT SES, SIS, SLS, SGHS, SHS, SEL, SIL, SLL, SGHL, SHL;
    DT dt;
    DT xnoh;    /* sgp4fix */
    dt = sat->t;
    /* sun terms */
    ZM=sat->deep.sun.zmo+g->nS*dt;
    ZF=ZM+2.*g->eS*SIN (ZM);
    SINZF=SIN (ZF);
    f2=.5*SINZF*SINZF-.25;
    f3=-.5*SINZF*COS (ZF);
    SES=sat->deep.sun.e[0]*f2+sat->deep.sun.e[1]*f3;
    SIS=sat->deep.sun.i[0]*f2+sat->deep.sun.i[1]*f3;
    SLS=sat->deep.sun.l[0]*f2+sat->deep.sun.l[1]*f3+sat->deep.sun.l[2]*SINZF;
    SGHS=sat->deep.sun.gh[0]*f2+sat->deep.sun.gh[1]*f3+sat->deep.sun.gh[2]*SINZF;
    SHS=sat->deep.sun.h[0]*f2+sat->deep.sun.h[1]*f3;
    /* moon terms */
    ZM=sat->deep.moon.zmo+g->nL*dt;
    ZF=ZM+2.*g->eL*SIN (ZM);
    SINZF=SIN (ZF);
    f2=.5*SINZF*SINZF-.25;
    f3=-.5*SINZF*COS (ZF);
    SEL=sat->deep.moon.e[0]*f2+sat->deep.moon.e[1]*f3;
    SIL=sat->deep.moon.i[0]*f2+sat->deep.moon.i[1]*f3;
    SLL=sat->deep.moon.l[0]*f2+sat->deep.moon.l[1]*f3+sat->deep.moon.l[2]*SINZF;
    SGHL=sat->deep.moon.gh[0]*f2+sat->deep.moon.gh[1]*f3+sat->deep.moon.gh[2]*SINZF;
    SHL=sat->deep.moon.h[0]*f2+sat->deep.moon.h[1]*f3;
    PE=SES+SEL;
    PINC=SIS+SIL;
    PL=SLS+SLL;

    PGH=SGHS+SGHL;
    PH=SHS+SHL;
    *XINC = *XINC+PINC;
    *EM = *EM+PE;
    sinis = SIN(*XINC);
    cosis = COS(*XINC);
    /*if (sat->o.i<.2) doubious sgp4fix */
    if (*XINC>=.2) {
        /* apply periodics directly */
        PH /= sinis;
        PGH -= cosis*PH;
        *OMGASM += PGH;
        *XNODES += PH;
        *XLL += PL;
    } else {
        /* apply periodics with lyddane modification */
        SINOK=SIN(*XNODES);
        COSOK=COS(*XNODES);
        ALFDP=sinis*SINOK;
        BETDP=sinis*COSOK;
        DALF=PH*COSOK+PINC*cosis*SINOK;
        DBET=-PH*SINOK+PINC*cosis*COSOK;
        ALFDP=ALFDP+DALF;
        BETDP=BETDP+DBET;
        *XNODES=FMOD(*XNODES, 2.*PI);
        XLS = *XLL+*OMGASM+cosis*(*XNODES);
        DLS=PL+PGH-PINC*(*XNODES)*sinis;
        XLS=XLS+DLS;
        xnoh = *XNODES;
        *XNODES=ATAN2(ALFDP,BETDP);
        if (FABS(xnoh-*XNODES)>PI) {
            if (*XNODES<xnoh)
                *XNODES += 2.*PI;
            else
                *XNODES -= 2.*PI;
        }
        *XLL = *XLL+PL;
        *OMGASM = XLS-*XLL-cosis*(*XNODES);
    }
}

