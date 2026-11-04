#ifndef PLACE_H
#define PLACE_H

#include "vecops.h"
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  double twist_per_repeat;
} PlaceOpt;

int place_parallel_transport(const Vec3 *pts, int K,
                             const Vec3 *unit_rel, size_t nAtoms,
                             Vec3 axis_unit,
                             const PlaceOpt *opt,
                             Vec3 *out);

#ifdef __cplusplus
}
#endif
#endif
