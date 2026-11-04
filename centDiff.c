#include "centDiff.h"

Vec3 tangent_central(const Vec3 *pts, int n, int i){
  if (!pts || n < 2) return (Vec3){0,0,0};
  if (i <= 0)     return vec3_nrm(vec3_sub(pts[1], pts[0]));
  if (i >= n - 1) return vec3_nrm(vec3_sub(pts[n-1], pts[n-2]));
  return vec3_nrm(vec3_sub(pts[i+1], pts[i-1]));
}
