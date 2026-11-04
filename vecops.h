#ifndef VECOPS_H
#define VECOPS_H

#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct { double x,y,z; } Vec3;

/* ---- basic ops ---- */
static inline Vec3 vec3_add(Vec3 a, Vec3 b){ return (Vec3){a.x+b.x, a.y+b.y, a.z+b.z}; }
static inline Vec3 vec3_sub(Vec3 a, Vec3 b){ return (Vec3){a.x-b.x, a.y-b.y, a.z-b.z}; }
static inline Vec3 vec3_scl(Vec3 a, double s){ return (Vec3){a.x*s, a.y*s, a.z*s}; }

static inline double vec3_dot(Vec3 a, Vec3 b){ return a.x*b.x + a.y*b.y + a.z*b.z; }

static inline Vec3 vec3_crs(Vec3 a, Vec3 b){
  return (Vec3){
    a.y*b.z - a.z*b.y,
    a.z*b.x - a.x*b.z,
    a.x*b.y - a.y*b.x
  };
}

static inline double vec3_len(Vec3 a){ return sqrt(vec3_dot(a,a)); }

static inline double vec3_len2(Vec3 a){ return vec3_dot(a,a); }

static inline double vec3_dist(Vec3 a, Vec3 b){ return vec3_len(vec3_sub(a,b)); }

static inline double vec3_dist2(Vec3 a, Vec3 b){ return vec3_len2(vec3_sub(a,b)); }

static inline Vec3 vec3_nrm(Vec3 a){
  double L = vec3_len(a); return (L > 0.0) ? vec3_scl(a, 1.0/L) : (Vec3){0,0,0}; }

/* ---- 3x3 matrix helpers ---- */
static inline void mat3_ident(double R[3][3]){
  R[0][0]=1; R[0][1]=0; R[0][2]=0;
  R[1][0]=0; R[1][1]=1; R[1][2]=0;
  R[2][0]=0; R[2][1]=0; R[2][2]=1;
}

static inline Vec3 mat3_apply(const double R[3][3], Vec3 v){
  return (Vec3){
    R[0][0]*v.x + R[0][1]*v.y + R[0][2]*v.z,
    R[1][0]*v.x + R[1][1]*v.y + R[1][2]*v.z,
    R[2][0]*v.x + R[2][1]*v.y + R[2][2]*v.z
  };
}

static inline void mat3_mul(const double A[3][3], const double B[3][3], double C[3][3]){
  for (int i=0;i<3;i++){
    for (int j=0;j<3;j++){
      C[i][j] = 0.0;
      for (int k=0;k<3;k++) C[i][j] += A[i][k]*B[k][j];
    }
  }
}

/* ---- rotations ---- */
static inline void rot_about_axis(Vec3 axis, double ang, double R[3][3]){
  axis = vec3_nrm(axis);
  double x=axis.x, y=axis.y, z=axis.z;
  double c = cos(ang), s = sin(ang), C = 1.0 - c;

  R[0][0] = c + x*x*C; R[0][1] = x*y*C - z*s; R[0][2] = x*z*C + y*s;
  R[1][0] = y*x*C + z*s; R[1][1] = c + y*y*C; R[1][2] = y*z*C - x*s;
  R[2][0] = z*x*C - y*s; R[2][1] = z*y*C + x*s; R[2][2] = c + z*z*C;
}

/* rotation matrix that maps v_from onto v_to */
static inline void rot_from_to(Vec3 v_from, Vec3 v_to, double R[3][3]){
  v_from = vec3_nrm(v_from);
  v_to   = vec3_nrm(v_to);

  double c = vec3_dot(v_from, v_to);
  Vec3 ax  = vec3_crs(v_from, v_to);
  double s = vec3_len(ax);

  if (s < 1e-12){
    if (c > 0.0){ mat3_ident(R); return; }
    Vec3 ortho = (fabs(v_from.x) < 0.9) ? (Vec3){1,0,0} : (Vec3){0,1,0};
    ax = vec3_nrm(vec3_crs(v_from, ortho));
    double x=ax.x, y=ax.y, z=ax.z;
    R[0][0] = -1 + 2*x*x; R[0][1] = 2*x*y;     R[0][2] = 2*x*z;
    R[1][0] = 2*y*x;     R[1][1] = -1 + 2*y*y; R[1][2] = 2*y*z;
    R[2][0] = 2*z*x;     R[2][1] = 2*z*y;     R[2][2] = -1 + 2*z*z;
    return;
  }

  ax = vec3_scl(ax, 1.0/s);
  double x=ax.x, y=ax.y, z=ax.z;
  double C = 1.0 - c;

  R[0][0] = c + x*x*C; R[0][1] = x*y*C - z*s; R[0][2] = x*z*C + y*s;
  R[1][0] = y*x*C + z*s; R[1][1] = c + y*y*C; R[1][2] = y*z*C - x*s;
  R[2][0] = z*x*C - y*s; R[2][1] = z*y*C + x*s; R[2][2] = c + z*z*C;
}

/*Special Primitive of central difference tanger for later use in curve fitting*/

/*static Vec3 tangent_central(const Vec3 *pts, int n, int i){
  if (n < 2) return (Vec3){0,0,0};
  if (i <= 0)     return vec3_nrm(vec3_sub(pts[1], pts[0]));
  if (i >= n - 1) return vec3_nrm(vec3_sub(pts[n-1], pts[n-2]));
  return vec3_nrm(vec3_sub(pts[i+1], pts[i-1]));
}
*/

static inline Vec3 vec3_proj_perp(Vec3 v, Vec3 axis_unit){
  /* remove component along axis_unit */
  return vec3_sub(v, vec3_scl(axis_unit, vec3_dot(v, axis_unit)));
}

/* Signed angle to rotate from_u -> to_u about axis_unit (right-hand rule).
   Requires from_u and to_u be perpendicular to axis_unit (or close). */
static inline double signed_angle_about_axis(Vec3 from_u, Vec3 to_u, Vec3 axis_unit){
  Vec3 c = vec3_crs(from_u, to_u);
  double s = vec3_dot(c, axis_unit);
  double d = vec3_dot(from_u, to_u);
  return atan2(s, d);
}

static inline Vec3 curve_discrete_normal(const Vec3 *pts, int n, int i){
	/* Returns a vector intended to point toward local center of curvature via second difference */
	if(!pts || n < 3) return (Vec3){0,0,0};
	if( i <= 0 || i >= n-1) return (Vec3){0,0,0};
	Vec3 a = vec3_sub(pts[i+1],pts[i]);
	Vec3 b = vec3_sub(pts[i],pts[i-1]);
	return vec3_sub(a,b);
}

static inline Vec3 discrete_dd(const Vec3 *p, int n, int i){
  if (!p || n < 3) return (Vec3){0,0,0};
  if (i <= 0 || i >= n-1) return (Vec3){0,0,0};
  Vec3 a = vec3_sub(p[i+1], p[i]);
  Vec3 b = vec3_sub(p[i], p[i-1]);
  return vec3_sub(a, b);
}

#ifdef __cplusplus
}
#endif
#endif
