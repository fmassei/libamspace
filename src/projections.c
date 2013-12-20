#include "projections.h"
#include "angles.h"

void amsp_projections_sph2latlon(DT out[2], const DT r[3])
{
    out[0] = r[0]*CRG;
    out[1] = r[1]*CRG;
}
void amsp_projections_latlon2mercator(DT out[2], const DT r[2], int w, int h)
{
    DT mn;
    out[0] = (r[1]+180.)*w/360.;
    mn = LOG(TAN((PI/4.)+(r[0]*CGR)/2.));
    out[1] = (h/2.)-(w*mn/(2.*PI));
}

void amsp_projections_latlon2utm(DT out[4], const DT r[3], gravconst_st *g)
{
    DT lat_r, lat, lon;
    DT coslat, sinlat, k0, a, b, e, n, A0, B0, C0, D0, E0, S, e1sq, nu, p,
       Ki, Kii, Kiii, Kiv, Kv;
    int zone, zone_cm;
    lat_r = (lat=r[0])*CGR;
    lon = r[1];
    coslat = COS(lat_r);
    sinlat = SIN(lat_r);
    k0 = 0.9996;
    a = g->rE*1000;
    b = g->prE*1000;
    e = SQRT(1.-P2(b/a));
    e1sq = e*e/(1-P2(e));
    n = (a-b)/(a+b);
    A0 = a*(1-n+(5.*n*n/4.)*(1.-n)+(81*P4(n)/64.)*(1.-n));
    B0 = (3*a*n/2.)*(1-n-(7*n*n/8.)*(1-n)+55*P4(n)/64.);
    C0 = (15*a*n*n/16.)*(1-n+(3*n*n/4.)*(1-n));
    D0 = (35*a*P3(n)/48.)*(1-n+11*n*n/16.);
    E0 = (315*a*P4(n)/51.)*(1-n);
    zone = 31+FLOOR(lon/6.);
    zone_cm = 6.*zone-183;
    p = (lon-zone_cm)*CGR;
    nu = a/SQRT(1-P2(e*sinlat));
    S = A0*lat_r - B0*sin(2.*lat_r) + C0*sin(4.*lat_r) -
        D0*sin(6.*lat_r) + E0*sin(8.*lat_r);
    Ki = S*k0;
    Kii = nu*sinlat*coslat/2.;
    Kiii = ((nu*sinlat*P3(coslat))/24.)*(5-P2(TAN(lat_r))+9.*e1sq*P2(coslat)+
            4.*P2(e1sq)*P4(coslat))*k0;
    Kiv = nu*coslat*k0;
    Kv = P3(coslat)*(nu/6.)*(1.-P2(TAN(lat_r))+e1sq*P2(coslat))*k0;
    out[0] = Ki+Kii*P2(p)+Kiii*P4(p);
    if (out[0]<0.) {
        out[3] = 0;
        out[0]+=10000000.;
    } else out[3] = 1;
    out[1] = 500000.+(Kiv*p+Kv*P3(p));
    out[2] = zone;
}

/*void amsp_projections_utm2map(DT out[2], const DT r[4], int w, int h)
{
    TODO
}*/

