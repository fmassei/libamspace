#ifndef H_PROJECTIONS_H
#define H_PROJECTIONS_H

#include "exports.h"
#include "precision.h"
#include "gravconst.h"

/* lat/lon */
LIB_PUBLIC void amsp_projections_sph2latlon(DT out[2], const DT r[3]);
LIB_PUBLIC void amsp_projections_latlon2mercator(DT out[2],
                                                 const DT r[2], int w, int h);
/* N/E/zone/isNorth */
LIB_PUBLIC void amsp_projections_latlon2utm(DT out[4], const DT r[3],
                                            gravconst_st *g);

#endif /* H_PROJECTIONS_H */
