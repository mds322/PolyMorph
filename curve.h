//curve Header
//
#ifndef CURVE_H
#define CURVE_H

#include "vecops.h"  /* Vec3 */
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef Vec3 (*CurvePosFn)(double t, void *user);

/* Workspace for parametric curve sampling (caller-owned memory). */
typedef struct {
  Vec3   *p;   /* length M */
  double *s;   /* length M, cumulative arc-length */
  int     M;   /* >=2 */
} CurveWS;

/* Build oversampled polyline + cumulative arc-length table for r(t), t in [t0,t1]. */
int curve_ws_build(CurvePosFn f, void *user, double t0, double t1, CurveWS *ws);

/* Resample from ws to N points equally spaced by arc-length. */
int curve_ws_resample_equal_arclen(const CurveWS *ws, int N, Vec3 *pts_out);

/* Resample from ws to ~step spacing; writes up to maxN, returns count or -1. */
int curve_ws_resample_by_step(const CurveWS *ws, double step, Vec3 *pts_out, int maxN);

/* ---- Parabola wrapper: y = a*x^2 + b*x + c, z=0.
   Parameter t is x itself (domain [t0,t1] == [x_min,x_max]). */
typedef struct {
  double a, b, c;
} Parabola;

Vec3 curve_pos_parabola(double t, void *user); /* user must point to Parabola */

/* ---- Utility: scaling sampled points ----
   Compute uniform scale factor so that the (median) distance between consecutive
   points is close to target_step (e.g., O4-O4 distance). */
double curve_scale_for_target_step(const Vec3 *pts, int N, double target_step);

/* Multiply all points by s (uniform XYZ scaling). */
void curve_scale_pts(Vec3 *pts, int N, double s);

/* ----Utility: Auto-Tune Scaling to avoid artifacts from over/undersampling range ---- */
/* Compute total arc-length over [t0,t1] by building ws (ws->p/ws->s must be allocated).
   Returns arc-length on success (>=0). Returns negative on error. */
double curve_ws_arclen(CurvePosFn f, void *user, double t0, double t1, CurveWS *ws);

/* Auto-tune t1 so that arc-length over [t0,t1] ~= target_len.
   Assumes t1_initial != t0. Uses expansion + bisection.
   On success writes *t1_io and returns 0. Returns -1 on failure. */
int curve_autotune_t1_for_target_len(CurvePosFn f, void *user,
                                     double t0, double *t1_io,
                                     double target_len,
                                     CurveWS *ws,
                                     int max_expand, int max_bisect,
                                     double rel_tol);

#ifdef __cplusplus
}
#endif
#endif
