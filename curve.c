//### `curve.c` (generic sampler + parabola position function)
//

#include "vecops.h"
#include "curve.h"
#include <math.h>
#include <stdlib.h> /* qsort */


static Vec3 v_lerp(Vec3 a, Vec3 b, double u){
  return (Vec3){ a.x + u*(b.x-a.x), a.y + u*(b.y-a.y), a.z + u*(b.z-a.z) };
}

static int lower_seg(const double *s, int M, double st){
  int lo=0, hi=M-1;
  while (hi - lo > 1){
    int mid = (lo+hi)/2;
    if (s[mid] <= st) lo = mid;
    else hi = mid;
  }
  return lo;
}

int curve_ws_build(CurvePosFn f, void *user, double t0, double t1, CurveWS *ws){
  if (!f || !ws || !ws->p || !ws->s) return -1;
  if (ws->M < 2) return -1;
  if (t1 == t0) return -1;

  int M = ws->M;
  double dt = (t1 - t0) / (M - 1);

  for (int i=0;i<M;i++){
    double t = t0 + i*dt;
    ws->p[i] = f(t, user);
  }

  ws->s[0] = 0.0;
  for (int i=1;i<M;i++){
    ws->s[i] = ws->s[i-1] + vec3_dist(ws->p[i], ws->p[i-1]);
  }
  return 0;
}

int curve_ws_resample_equal_arclen(const CurveWS *ws, int N, Vec3 *pts_out){
  if (!ws || !ws->p || !ws->s || !pts_out) return -1;
  if (ws->M < 2 || N < 2) return -1;

  int M = ws->M;
  double L = ws->s[M-1];

  if (!(L > 0.0)){
    for (int i=0;i<N;i++) pts_out[i] = ws->p[0];
    return 0;
  }

  double ds = L / (N - 1);

  pts_out[0]   = ws->p[0];
  pts_out[N-1] = ws->p[M-1];

  for (int i=1;i<N-1;i++){
    double st = i * ds;
    int k = lower_seg(ws->s, M, st);
    double denom = ws->s[k+1] - ws->s[k];
    double u = (denom > 0.0) ? (st - ws->s[k]) / denom : 0.0;
    pts_out[i] = v_lerp(ws->p[k], ws->p[k+1], u);
  }
  return 0;
}

int curve_ws_resample_by_step(const CurveWS *ws, double step, Vec3 *pts_out, int maxN){
  if (!ws || !ws->p || !ws->s || !pts_out) return -1;
  if (ws->M < 2 || maxN < 2) return -1;
  if (!(step > 0.0)) return -1;

  int M = ws->M;
  double L = ws->s[M-1];

  if (!(L > 0.0)){
    pts_out[0] = ws->p[0];
    pts_out[1] = ws->p[M-1];
    return 2;
  }

  int n = 0;
  double st = 0.0;

  while (st <= L + 1e-12 && n < maxN){
    int k = lower_seg(ws->s, M, st);
    double denom = ws->s[k+1] - ws->s[k];
    double u = (denom > 0.0) ? (st - ws->s[k]) / denom : 0.0;
    pts_out[n++] = v_lerp(ws->p[k], ws->p[k+1], u);
    st += step;
  }

  if (n < maxN){
    pts_out[n++] = ws->p[M-1];
  } else {
    pts_out[maxN-1] = ws->p[M-1];
    n = maxN;
  }
  return n;
}

/* ---- Parabola wrapper ----
   Parameter t is x. Domain [t0,t1] passed to curve_ws_build is [x_min,x_max]. */
Vec3 curve_pos_parabola(double t, void *user){
  const Parabola *P = (const Parabola*)user;
  double x = t;
  double y = P->a*x*x + P->b*x + P->c;
  return (Vec3){ x, y, 0.0 };
}

/* ---- scaling helpers ---- */

static int cmp_dbl(const void *a, const void *b){
  double da = *(const double*)a, db = *(const double*)b;
  return (da < db) ? -1 : (da > db) ? 1 : 0;
}

double curve_scale_for_target_step(const Vec3 *pts, int N, double target_step){
  if (!pts || N < 2) return 1.0;
  if (!(target_step > 0.0)) return 1.0;

  int m = N - 1;

  enum { MAX_M = 4096 };
  double d[MAX_M];

  int stride = 1;
  int use_m = m;
  if (use_m > MAX_M){
    stride = (use_m + MAX_M - 1) / MAX_M; /* ceil */
    use_m = (m + stride - 1) / stride;
    if (use_m > MAX_M) use_m = MAX_M;
  }

  int j = 0;
  for (int i=0; i<m && j<use_m; i += stride){
    d[j++] = vec3_dist(pts[i+1], pts[i]);
  }
  use_m = j;
  if (use_m < 1) return 1.0;

  qsort(d, (size_t)use_m, sizeof(d[0]), cmp_dbl);

  double med = (use_m % 2) ? d[use_m/2] : 0.5*(d[use_m/2 - 1] + d[use_m/2]);
  if (!(med > 0.0)) return 1.0;

  return target_step / med;
}

void curve_scale_pts(Vec3 *pts, int N, double s){
  if (!pts) return;
  for (int i=0;i<N;i++){
    pts[i].x *= s;
    pts[i].y *= s;
    pts[i].z *= s;
  }
}
double curve_ws_arclen(CurvePosFn f, void *user, double t0, double t1, CurveWS *ws){
  if (curve_ws_build(f, user, t0, t1, ws) != 0) return -1.0;
  if (!ws || !ws->s || ws->M < 2) return -1.0;
  return ws->s[ws->M - 1];
}

int curve_autotune_t1_for_target_len(CurvePosFn f, void *user,
                                     double t0, double *t1_io,
                                     double target_len,
                                     CurveWS *ws,
                                     int max_expand, int max_bisect,
                                     double rel_tol)
{
  if (!t1_io || !ws) return -1;
  if (!(target_len > 0.0)) return -1;
  if (!(rel_tol > 0.0)) rel_tol = 1e-3;
  if (max_expand <= 0) max_expand = 32;
  if (max_bisect <= 0) max_bisect = 64;

  double t1 = *t1_io;
  if (t1 == t0) return -1;

  /* direction */
  double dir = (t1 > t0) ? 1.0 : -1.0;

  /* current length */
  double L = curve_ws_arclen(f, user, t0, t1, ws);
  if (!(L >= 0.0)) return -1;

  double Ltarget = target_len;
  double tol = rel_tol * Ltarget;

  /* If already close, keep */
  if (fabs(L - Ltarget) <= tol){
    *t1_io = t1;
    return 0;
  }

  /* Bracket solution in t1: find [lo, hi] such that L(lo) <= Ltarget <= L(hi) */
  double lo = t0, hi = t1;
  double Llo = 0.0, Lhi = L;

  if (L < Ltarget){
    /* need to expand away from t0 */
    lo = t1;
    Llo = L;
    hi = t1;
    Lhi = L;

    double step = (t1 - t0);
    if (step == 0.0) step = dir * 1.0;

    for (int it=0; it<max_expand; it++){
      hi = hi + step;        /* extend */
      step *= 2.0;

      Lhi = curve_ws_arclen(f, user, t0, hi, ws);
      if (!(Lhi >= 0.0)) return -1;
      if (Lhi >= Ltarget) break;
    }
    if (Lhi < Ltarget) return -1;
    /* now bracket is [lo=t1_initial, hi] with Llo < Ltarget <= Lhi */
    lo = *t1_io;
    Llo = L;
  } else {
    /* L > Ltarget: we need to move t1 toward t0 to reduce length.
       We'll search between t0 and current t1. */
    lo = t0;
    Llo = 0.0;
    hi = t1;
    Lhi = L;

    /* If even a very small interval overshoots (shouldn’t happen unless Ltarget tiny),
       bisection will still handle it. */
  }

  /* Bisection on t1 */
  for (int it=0; it<max_bisect; it++){
    double mid = 0.5*(lo + hi);
    if (mid == lo || mid == hi) break;

    double Lmid = curve_ws_arclen(f, user, t0, mid, ws);
    if (!(Lmid >= 0.0)) return -1;

    double err = Lmid - Ltarget;
    if (fabs(err) <= tol){
      *t1_io = mid;
      return 0;
    }

    if (err < 0.0){
      lo = mid;
      Llo = Lmid;
    } else {
      hi = mid;
      Lhi = Lmid;
    }
  }

  /* best effort: choose closer of lo/hi */
  double err_lo = fabs(Llo - Ltarget);
  double err_hi = fabs(Lhi - Ltarget);
  *t1_io = (err_lo <= err_hi) ? lo : hi;
  return 0;
}
