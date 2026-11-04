#include "vecops.h"
#include "curvature.h"
#include <math.h>

double polyline_min_circumradius(const Vec3 *p, int N){
  if (!p || N < 3) return HUGE_VAL;

  double Rmin = HUGE_VAL;

  for (int i=1; i<N-1; i++){
    Vec3 a = vec3_sub(p[i],   p[i-1]);
    Vec3 b = vec3_sub(p[i+1], p[i]);
    Vec3 c = vec3_sub(p[i+1], p[i-1]);

    double la2 = vec3_len2(a);
    double lb2 = vec3_len2(b);
    double lc2 = vec3_len2(c);
    if (!(la2 > 0.0) || !(lb2 > 0.0) || !(lc2 > 0.0)) continue;

    double la = sqrt(la2), lb = sqrt(lb2), lc = sqrt(lc2);

    /* area2 = |a x b| = 2*triangle_area */
    Vec3 axb = vec3_crs(a, b);
    double area2 = vec3_len(axb);
    if (!(area2 > 1e-15)) continue; /* collinear */

    double A = 0.5 * area2;
    double R = (la * lb * lc) / (4.0 * A);
    if (R < Rmin) Rmin = R;
  }

  return Rmin;
}
