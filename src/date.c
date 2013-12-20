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
#include "date.h"
#include "angles.h"
#include <math.h>

DT amsp_date_tm2julian(const dtm_st *dtm)
{
    return 367.*(dtm->tm.tm_year+1900)-
       FLOOR((7*(dtm->tm.tm_year+1900+FLOOR((dtm->tm.tm_mon+1+9)/12.)))*.25)+
       FLOOR(275*(dtm->tm.tm_mon+1)/9.)+
       dtm->tm.tm_mday+1721013.5+
       ((dtm->sec/60.+dtm->tm.tm_min)/60.+dtm->tm.tm_hour)/24.;
}

void amsp_date_YD2tm(dtm_st *out, int year, DT days)
{
    int i, inttmp, yday;
    double tmp;
    int dm[] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
    out->tm.tm_year = year-1900;
    yday = (int)FLOOR(days);
    if (year%4==0) dm[1] = 29;
    i = 0;
    inttmp = 0;
    while ((yday>inttmp+dm[i])&&(i<11)) {
        inttmp += dm[i];
        ++i;
    }
    out->tm.tm_mon = i;
    out->tm.tm_mday = yday-inttmp;

    tmp = (days-yday)*24.;
    out->tm.tm_hour = (int)FLOOR(tmp);
    tmp = (tmp-out->tm.tm_hour)*60.;
    out->tm.tm_min = (int)FLOOR(tmp);
    out->sec = (tmp-out->tm.tm_min)*60.;
}

void amsp_date_julian2tm(dtm_st *out, DT j)
{
    int ly;
    DT d, t, tmp;
    tmp = j-2415019.5;
    t = tmp/365.25;
    out->tm.tm_year = (int)FLOOR(t);
    ly = (int)FLOOR((out->tm.tm_year-1)*.25);

    d = tmp - (out->tm.tm_year*365.0+ly)+0.00000000001;
    if (d<1.) {
        --out->tm.tm_year;
        ly = (int)FLOOR((out->tm.tm_year-1)*.25);
        d = tmp - ((out->tm.tm_year)*365.0+ly);
    }
    amsp_date_YD2tm(out, out->tm.tm_year+1900, d);
    out->sec -= 0.00000086400;
}

DT amsp_date_gstime(DT j)
{
    DT tmp, tut1;
    tut1 = (j-2451545.)/36525.;
    tmp = 67310.54841 + (876600.0*3600.+8640184.812866)*tut1 +
         0.093104*P2(tut1) - 6.2e-6*P3(tut1);
    tmp = FMOD(tmp*CGR/240., 2.*PI);
    if (tmp<0.)
        tmp += 2.*PI;
    return tmp;
}

/*#include <stdio.h>
int main(void)
{
    struct tm tm;
    time_t t;
    t = time(NULL);
    localtime_r(&t, &tm);
    printf("%lf\n", toJulian(&tm));
    return 0;
}*/
