#ifndef CURVATURE_H
#define CURVATURE_H

#include "vecops.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Minimum circumradius over polyline vertices i=1..N-2.
   Returns HUGE_VAL if N<3 or all vertices are collinear/degenerate. */
double polyline_min_circumradius(const Vec3 *p, int N);

#ifdef __cplusplus
}
#endif
#endif
