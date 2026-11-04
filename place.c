#include "centDiff.h"
#include "place.h"
#include <string.h>
#include <math.h>

int place_parallel_transport(const Vec3 *pts, int K,
                             const Vec3 *unit_rel, size_t nAtoms,
                             Vec3 axis_unit,
                             const PlaceOpt *opt,
                             Vec3 *out)
{
  if (!pts || !unit_rel || !out) return -1;
  if (K < 1 || nAtoms == 0) return -1;

  double twist_per_repeat = 0.0;
  if (opt) twist_per_repeat = opt->twist_per_repeat;

  Vec3 T_prev = tangent_central(pts, K, 0);
  if (vec3_len(T_prev) < 1e-12) T_prev = axis_unit;

  double Q[3][3];
  rot_from_to(axis_unit, T_prev, Q);

  Vec3 up_prev;
  {
    Vec3 dd0 = discrete_dd(pts, K, (K > 2) ? 1 : 0);
    dd0 = vec3_proj_perp(dd0, T_prev);
    if (vec3_len(dd0) < 1e-12){
      Vec3 ortho = (fabs(T_prev.x) < 0.9) ? (Vec3){1,0,0} : (Vec3){0,1,0};
      dd0 = vec3_crs(T_prev, ortho);
    }
    up_prev = vec3_nrm(dd0);
  }

  for (int r=0; r<K; r++){
    Vec3 T = tangent_central(pts, K, r);
    if (vec3_len(T) < 1e-12) T = T_prev;

    if (r > 0){
      double R_step[3][3], Q_new[3][3];
      rot_from_to(T_prev, T, R_step);
      mat3_mul(R_step, Q, Q_new);
      memcpy(Q, Q_new, sizeof(Q));

      up_prev = mat3_apply(R_step, up_prev);
      up_prev = vec3_proj_perp(up_prev, T);
      if (vec3_len(up_prev) < 1e-12){
        Vec3 ortho = (fabs(T.x) < 0.9) ? (Vec3){1,0,0} : (Vec3){0,1,0};
        up_prev = vec3_crs(T, ortho);
      }
      up_prev = vec3_nrm(up_prev);
    }

    Vec3 up_des = discrete_dd(pts, K, r);
    up_des = vec3_proj_perp(up_des, T);

    if (vec3_len(up_des) < 1e-8){
      up_des = up_prev;
    } else {
      up_des = vec3_nrm(up_des);
      up_prev = up_des;
    }

    Vec3 up_cur = mat3_apply(Q, (Vec3){1,0,0});
    up_cur = vec3_proj_perp(up_cur, T);
    if (vec3_len(up_cur) < 1e-12){
      up_cur = mat3_apply(Q, (Vec3){0,1,0});
      up_cur = vec3_proj_perp(up_cur, T);
    }
    if (vec3_len(up_cur) < 1e-12){
      up_cur = up_des;
    } else {
      up_cur = vec3_nrm(up_cur);
    }

    double theta_lock = signed_angle_about_axis(up_cur, up_des, T);
    double theta = theta_lock + r * twist_per_repeat;

    double R_twist[3][3], R_total[3][3];
    rot_about_axis(T, theta, R_twist);
    mat3_mul(R_twist, Q, R_total);

    for (size_t i=0;i<nAtoms;i++){
      Vec3 pr = mat3_apply(R_total, unit_rel[i]);
      out[(size_t)r*nAtoms + i] = vec3_add(pr, pts[r]);
    }

    T_prev = T;
  }

  return 0;
}
