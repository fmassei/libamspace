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
#ifndef H_DEEP_H
#define H_DEEP_H

#include "exports.h"
#include "precision.h"
#include "sgpcommon.h"
#include "gravconst.h"

/* as a reference: variable names in the "restriced four body..." paper:
 * w = argp, omega = node, gamma = ma, n = mm, i = incl,
 * G = argp, H = node, N = mm
 */

/* body components, for moon and sun, used by deep space functions.
 * none has a specifical meaning, just intermediate values for making spherical
 * stuff a little bit easier. */
typedef struct body_s body_st;
struct body_s {
    DT e[2], i[2], l[3], gh[3], h[2], zmo;
};

/* type of orbital resonance */
typedef enum resonanceType_e {
    RESONANCE_NONE,
    RESONANCE_ONEDAY,
    RESONANCE_HALFDAY,
} resonanceType_et;

/* deep space stuff */
typedef struct dspace_s dspace_st;
struct dspace_s {
    resonanceType_et res;
    DT q[6][5][4][4];   /* (Q_lmpq) coefficients for the four-body integration
                            (amplitude of term in spherical harmonic potential
                            disturbing function) where l=degree of spherical
                            harmonic, m=order of spherical harmonic, p=index
                            of inclination function, q=index of eccentricity
                            function */
    DT gst;             /* epoch GST */
    DT el[3];           /* ? */
    DT lam, fact;       /* start/step for the integrator */
    body_st moon, sun;  /* body tmp variables */
    orbitalelems_st de; /* sun/moon effects to apply on derivatives */
};

LIB_LOCAL void dpinit(sat_st *sat, const gravconst_st *g);
LIB_LOCAL void dpsec(sat_st *sat,
            DT *XLL, DT *OMGASM, DT *XNODES, DT *EM, DT *XINC, DT *XN);
LIB_LOCAL void dpper(sat_st *sat, const gravconst_st *g,
            DT *EM, DT *XINC, DT *OMGASM, DT *XNODES, DT *XLL);

#endif /* H_DEEP_H */
