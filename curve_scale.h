#ifndef CURVE_SCALE_H
#define CURVE_SCALE_H

#include "vecops.h"
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  double scale_step;
  double scale_curv;
  double scale_chosen;
  double r_unit;
  double Rmin;
} CurveScaleReport;

double unit_bounding_radius(const Vec3 *unit_rel, size_t nAtoms);

double curve_choose_global_scale(const Vec3 *pts, int K,
                                 const Vec3 *unit_rel, size_t nAtoms,
                                 double d_target,
                                 double curv_factor,
                                 CurveScaleReport *rep);

#ifdef __cplusplus
}
#endif
#endif
