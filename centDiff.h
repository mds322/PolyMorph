#ifndef CENTDIFF_H
#define CENTDIFF_H

#include "vecops.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Central-difference tangent for a polyline. Returns unit vector or {0,0,0}. */
Vec3 tangent_central(const Vec3 *pts, int n, int i);

#ifdef __cplusplus
}
#endif
#endif
