#include "curve_scale.h"
#include "curve.h"       /* curve_scale_for_target_step */
#include "curvature.h"   /* polyline_min_circumradius */
#include <math.h>

double unit_bounding_radius(const Vec3 *unit_rel, size_t nAtoms){
  if (!unit_rel || nAtoms == 0) return 0.0;
  double r = 0.0;
  for (size_t i=0;i<nAtoms;i++){
    double a = vec3_len(unit_rel[i]);
    if (a > r) r = a;
  }
  return r;
}

double curve_choose_global_scale(const Vec3 *pts, int K,
                                 const Vec3 *unit_rel, size_t nAtoms,
                                 double d_target,
                                 double curv_factor,
                                 CurveScaleReport *rep)
{
  CurveScaleReport tmp;
  tmp.scale_step = 1.0;
  tmp.scale_curv = 1.0;      /* NOTE: now means "extra factor beyond step" */
  tmp.scale_chosen = 1.0;    /* total factor = scale_step * scale_curv */
  tmp.r_unit = unit_bounding_radius(unit_rel, nAtoms);
  if (!(tmp.r_unit > 0.0)) tmp.r_unit = 1.0;

  /* curvature of the *unscaled* resampled polyline */
  tmp.Rmin = polyline_min_circumradius(pts, K);

  /* step-based scale: bring median step to d_target (in Angstrom) */
  tmp.scale_step = curve_scale_for_target_step(pts, K, d_target);
  if (!(tmp.scale_step > 0.0)) tmp.scale_step = 1.0;

  if (!(curv_factor > 0.0)) curv_factor = 3.0;

  /* curvature-based scale must be evaluated in final length units:
     circumradius scales linearly with uniform scaling. */
  tmp.scale_curv = 1.0;
  if (isfinite(tmp.Rmin) && tmp.Rmin > 0.0 && tmp.Rmin < HUGE_VAL/2){
    double R_target = curv_factor * tmp.r_unit;
    double R_after_step = tmp.Rmin * tmp.scale_step;
    if (R_after_step < R_target){
      /* additional factor beyond step scaling */
      tmp.scale_curv = R_target / R_after_step;
    }
  }

  /* total scale */
  tmp.scale_chosen = tmp.scale_step * tmp.scale_curv;

  if (rep) *rep = tmp;
  return tmp.scale_chosen;
}
