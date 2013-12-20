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
#include "elp82.h"
#include "parseutils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct elp82data_s {
    DT *pc[3][6];
    DT *per[3][11][3];
    DT w[3][5];
    int nterm[3][12];
};

static int countLines(FILE *f)
{
    int lines = 1;
    int c;
    while ((c = getc(f))!=EOF)
        if (c=='\n') ++lines;
    rewind(f);
    return lines;
}

void amsp_elp82_destroyData(elp82data_st **data)
{
    int i, j, k;
    if (data==NULL || *data==NULL) return;
    for (i=0; i<3; ++i) {
        for (j=0; j<6; ++j)
            free((*data)->pc[i][j]);
        for (j=0; j<11; ++j)
            for (k=0; k<3; ++k)
                free((*data)->per[i][j][k]);
    }
    free(*data);
    data = NULL;
}

elp82data_st *amsp_elp82_getData(const char *datadir)
{
    char *fname = NULL;
    char buf[150];
    FILE *in = NULL;
    int addSlash = (datadir[strlen(datadir)-1]!='/');
    int ifile;
    elp82data_st *data;

    DT cpi = PI,
       pis2 = cpi/2.,
       rad = 648000./cpi,
       deg = cpi/180.,
       c1 = 60.,
       c2 = 3600.;
    DT am = 0.074801329518,
       dtasm = 2.*0.002571881335/(3.*am);
    DT eart[5], peri[5], p[8][2], del[4][5], zeta[2], precess;
    DT delnu, delg, delnp, dele, delep;
    int di;

    if ((data = calloc(1, sizeof(*data)))==NULL)
        return NULL;

    data->w[0][0] = (218.+18./c1+59.95571/c2)*deg;
    data->w[1][0] = (83.+21./c1+11.67475/c2)*deg;
    data->w[2][0] = (125.+2./c1+40.39816/c2)*deg;
    eart[0] = (100.+27./c1+59.22059/c2)*deg;
    peri[0] = (102.+56./c1+14.42753/c2)*deg;
    data->w[0][1] = 1732559343.73604/rad;
    data->w[1][1] = 14643420.2632/rad;
    data->w[2][1] = -6967919.3622/rad;
    eart[1] = 129597742.2758/rad;
    peri[1] = 1161.2283/rad;
    data->w[0][2] = -5.8883/rad;
    data->w[1][2] = -38.2776/rad;
    data->w[2][2] = 6.3622/rad;
    eart[2] = -0.0202/rad;
    peri[2] = 0.5327/rad;
    data->w[0][3] = 0.6604e-2/rad;
    data->w[1][3] = -0.45047e-1/rad;
    data->w[2][3] = 0.7625e-2/rad;
    eart[3] = 0.9e-5/rad;
    peri[3] = -0.138e-3/rad;
    data->w[0][4] = -0.3169e-4/rad;
    data->w[1][4] = 0.21301e-3/rad;
    data->w[2][4] = -0.3586e-4/rad;
    eart[4] = 0.15e-6/rad;
    peri[4] = 0.;

    precess=5029.0966/rad;

    p[0][0] = (252.+15./c1+3.25986/c2)*deg;
    p[1][0] = (181.+58./c1+47.28305/c2)*deg;
    p[2][0] = eart[0];
    p[3][0] = (355.+25./c1+59.78866/c2)*deg;
    p[4][0] = (34.+21./c1+5.34212/c2)*deg;
    p[5][0] = (50.+4./c1+38.89694/c2)*deg;
    p[6][0] = (314.+3./c1+18.01841/c2)*deg;
    p[7][0] = (304.+20./c1+55.19575/c2)*deg;
    p[0][1] = 538101628.68898/rad;
    p[1][1] = 210664136.43355/rad;
    p[2][1] = eart[1];
    p[3][1] = 68905077.59284/rad;
    p[4][1] = 10925660.42861/rad;
    p[5][1] = 4399609.65932/rad;
    p[6][1] = 1542481.19393/rad;
    p[7][1] = 786550.32074/rad;

    delnu=+0.55604/rad/data->w[0][1];
    dele=+0.01789/rad;
    delg=-0.08066/rad;
    delnp=-0.06424/rad/data->w[0][1];
    delep=-0.12879/rad;

    for (di=0; di<5; ++di) {
        del[0][di] = data->w[0][di]-eart[di];
        del[3][di] = data->w[0][di]-data->w[2][di];
        del[2][di] = data->w[0][di]-data->w[1][di];
        del[1][di] = eart[di]-peri[di];
    }
    del[0][0] = del[0][0]+cpi;
    zeta[0] = data->w[0][0];
    zeta[1] = data->w[0][1]+precess;
    
    for (ifile=0; ifile<36; ++ifile) {
        int ir, itab, iv, nrec;
        ir = 0;
        itab = ifile/3;
        iv = ifile%3;
        if ((fname = malloc(strlen(datadir)+addSlash+6))==NULL) {
            goto badexit;
        }
        sprintf(fname, "%s", datadir);
        if (addSlash) strcat(fname, "/");
        sprintf(fname, "%sELP%d", fname, ifile+1);
        if ((in = fopen(fname, "rb"))==NULL) {
            perror("fopen");
            goto badexit;
        }
        nrec = countLines(in)-1;
        if (nrec<=1) {
            goto badexit;
        }
        if (ifile<3) {
            int i;
            for (i=0; i<6; ++i)
                if ((data->pc[iv][i] = calloc(nrec, sizeof(DT)))==NULL)
                    goto badexit;
        } else {
            int i;
            for (i=0; i<3; ++i)
                if ((data->per[iv][itab-1][i] =
                            calloc(nrec, sizeof(DT)))==NULL)
                    goto badexit;
        }
        fgets(buf, sizeof(buf), in);    /* skip file description */
        while (fgets(buf, sizeof(buf), in)!=NULL) {
            if (buf[strlen(buf)-1]!='\n' && buf[strlen(buf)-1]!='\r')
                goto badexit;
            if (ifile<3) {
                int ilu[4];
                DT coef[7], zone[6], tgv, y;
                int i, k, sk = 0;
                for (i=0; i<4; ++i, sk+=3)
                    if (parseutils_buf2int(&ilu[i], buf+sk, 3)!=0)
                        goto badexit;
                if (parseutils_buf2double(&coef[0], buf+sk, 15)!=0)
                    goto badexit;
                sk += 15;
                for (i=1; i<7; ++i, sk+=12)
                    if (parseutils_buf2double(&coef[i], buf+sk, 12)!=0)
                        goto badexit;
                tgv = coef[1]+dtasm*coef[5];
                if (ifile==2) coef[0] = coef[0]-2.*coef[0]*delnu/3.;
                zone[0] = coef[0]+tgv*(delnp-am*delnu)+coef[2]*delg+
                          coef[3]*dele+coef[4]*delep;
                for (k=0; k<5; ++k) {
                    y = 0.;
                    for (i=0; i<4; ++i)
                        y += ilu[i]*del[i][k];
                    zone[k+1] = y;
                }
                if (iv==2) zone[1]+=pis2;
                for (i=0; i<6; ++i)
                    data->pc[iv][i][ir] = zone[i];
            } else if ((ifile>=3&&ifile<9)||(ifile>=21&&ifile<36)) {
                int ilu[5];
                DT pha, xx, zone[3], y;
                int i, k, sk = 0;
                for (i=0; i<5; ++i, sk+=3)
                    if (parseutils_buf2int(&ilu[i], buf+sk, 3)!=0)
                        goto badexit;
                if (parseutils_buf2double(&pha, buf+sk, 10)!=0)
                    goto badexit;
                sk += 10;
                if (parseutils_buf2double(&xx, buf+sk, 10)!=0)
                    goto badexit;
                zone[0] = xx;
                for (k=0; k<2; ++k) {
                    y = (k==0) ? (pha*deg) : 0.;
                    y += ilu[0]*zeta[k];
                    for (i=0; i<4; ++i)
                        y += ilu[i+1]*del[i][k];
                    zone[k+1] = y;
                }
                for (i=0; i<3; ++i)
                    data->per[iv][itab-1][i][ir] = zone[i];
            } else {
                int ipla[11];
                DT pha, xx, zone[3], y;
                int i, k, sk = 0;
                for (i=0; i<11; ++i, sk+=3)
                    if (parseutils_buf2int(&ipla[i], buf+sk, 3)!=0)
                        goto badexit;
                if (parseutils_buf2double(&pha, buf+sk, 10)!=0)
                    goto badexit;
                sk += 10;
                if (parseutils_buf2double(&xx, buf+sk, 10)!=0)
                    goto badexit;
                zone[0] = xx;
                if (ifile<15) {
                    for (k=0; k<2; ++k) {
                        y = (k==0) ? (pha*deg) : 0.;
                        y += ipla[8]*del[0][k]+ipla[9]*del[2][k]+
                             ipla[10]*del[3][k];
                        for (i=0; i<8; ++i)
                            y += ipla[i]*p[i][k];
                        zone[k+1] = y;
                    }
                } else {
                    for (k=0; k<2; ++k) {
                        y = (k==0) ? (pha*deg) : 0.;
                        for (i=0; i<4; ++i)
                            y += ipla[i+7]*del[i][k];
                        for (i=0; i<7; ++i)
                            y += ipla[i]*p[i][k];
                        zone[k+1] = y;
                    }
                }
                for (i=0; i<3; ++i)
                    data->per[iv][itab-1][i][ir] = zone[i];
            }
            ++ir;
        }
        data->nterm[iv][itab] = ir;
        fclose(in);
        free(fname);
        in = NULL;
        fname = NULL;
    }
    return data;
badexit:
    if (in) fclose(in);
    free(fname);
    amsp_elp82_destroyData(&data);
    return NULL;
}

void amsp_elp82_getPos(const elp82data_st *data, double jdate, double out[3])
{
    DT t[4], x0, x1, x2;
    DT p1, p2, p3, p4, p5, q1, q2, q3, q4, q5;
    DT x, y;
    int iv, itab, nt;
    DT pw, qw, ra, pwqw, pw2, qw2;
    DT ath = 384747.9806743165,
       a0 = 384747.9806448954,
       rad = 648000./PI;

    p1=0.10180391e-4;
    p2=0.47020439e-6;
    p3=-0.5417367e-9;
    p4=-0.2507948e-11;
    p5=0.463486e-14;
    q1=-0.113469002e-3;
    q2=0.12372674e-6;
    q3=0.1265417e-8;
    q4=-0.1371808e-11;
    q5=-0.320334e-14; 

    t[0] = (jdate-2451545.)/36525.;
    t[1] = t[0]*t[0]; t[2] = t[1]*t[0]; t[3] = t[2]*t[0];
    for (iv=0; iv<3; ++iv) {
        out[iv] = 0.;
        for (itab=0; itab<12; ++itab) {
            for (nt=0; nt<data->nterm[iv][itab]; ++nt) {
                if (itab==0) {
                    int k;
                    x = data->pc[iv][0][nt];
                    y = data->pc[iv][1][nt];
                    for (k=0; k<4; ++k)
                        y += data->pc[iv][k+2][nt]*t[k];
                } else {
                    x = data->per[iv][itab-1][0][nt];
                    y = data->per[iv][itab-1][1][nt] +
                        data->per[iv][itab-1][2][nt]*t[0];
                }
                if (itab==2 || itab==4 || itab==6 || itab==8)
                    x*=t[0];
                else if (itab==11)
                    x*=t[1];
                out[iv] += x*SIN(y);
            }
        }
    }
    out[0] = out[0]/rad + data->w[0][0] + 
             data->w[0][1]*t[0] + data->w[0][2]*t[1] + data->w[0][3]*t[2]+
             data->w[0][4]*t[3];
    out[1] /= rad;
    out[2] *= a0/ath;

    x0 = out[2]*COS(out[1]);
    x1 = x0*SIN(out[0]);
    x0 = x0*COS(out[0]);
    x2 = out[2]*SIN(out[1]);

    pw = (p1+p2*t[0]+p3*t[1]+p4*t[2]+p5*t[3])*t[0];
    qw = (q1+q2*t[0]+q3*t[1]+q4*t[2]+q5*t[3])*t[0];
    ra = 2.*SQRT(1.-P2(pw)-P2(qw));
    pwqw = 2.*pw*qw;
    pw2 = 1.-2*P2(pw);
    qw2 = 1.-2*P2(qw);
    pw *= ra;
    qw *= ra;

    out[0] = pw2*x0+pwqw*x1+pw*x2;
    out[1] = pwqw*x0+qw2*x1-qw*x2;
    out[2] = -pw*x0+qw*x1+(pw2+qw2-1.)*x2;
}

