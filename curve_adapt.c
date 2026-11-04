#include "curve_adapt.h"
#include <math.h>

static double clamp1(double x){
  if (x < -1.0) return -1.0;
  if (x >  1.0) return  1.0;
  return x;
}

static double vdot(Vec3 a, Vec3 b){ return a.x*b.x + a.y*b.y + a.z*b.z; }
static double vlen(Vec3 a){ return sqrt(vdot(a,a)); }
static Vec3 vsub(Vec3 a, Vec3 b){ return (Vec3){a.x-b.x, a.y-b.y, a.z-b.z}; }
static Vec3 vnrm(Vec3 a){
  double L = vlen(a);
  return (L > 0.0) ? (Vec3){a.x/L, a.y/L, a.z/L} : (Vec3){0,0,0};
}

double curve_ws_max_turn_angle(const CurveWS *ws){
  if (!ws || !ws->p || ws->M < 3) return -1.0;

  int M = ws->M;
  double th_max = 0.0;

  Vec3 v0 = vsub(ws->p[1], ws->p[0]);
  double L0 = vlen(v0);
  if (!(L0 > 0.0)) return 0.0;
  Vec3 t_prev = vnrm(v0);

  for (int i=1; i<M-1; i++){
    Vec3 v1 = vsub(ws->p[i+1], ws->p[i]);
    double L1 = vlen(v1);
    if (!(L1 > 0.0)) continue;

    Vec3 t = vnrm(v1);
    double c = clamp1(vdot(t_prev, t));
    double th = acos(c);
    if (th > th_max) th_max = th;
    t_prev = t;
  }

  return th_max;
}

int curve_ws_build_adaptive_turn(CurvePosFn f, void *user,
                                 double t0, double t1,
                                 CurveWS *ws,
                                 int maxM,
                                 double theta_limit_rad)
{
  if (!ws || !ws->p || !ws->s || ws->M < 2) return -1;
  if (!(theta_limit_rad > 0.0)) theta_limit_rad = (10.0 * M_PI / 180.0); /* default 10 deg */
  if (maxM < ws->M) maxM = ws->M;

  if (curve_ws_build(f, user, t0, t1, ws) != 0) return -1;

  if (ws->M < 3) return 0;

  double th_max = curve_ws_max_turn_angle(ws);
  if (!(th_max >= 0.0)) return -1;

  if (th_max <= theta_limit_rad) return 0;

  int newM = ws->M * 2;
  if (newM > maxM) newM = maxM;
  if (newM <= ws->M) return 0; /* can't improve further; accept */
  ws->M = newM;
  return +1; /* caller must realloc buffers and retry */
}
