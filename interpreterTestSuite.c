/* test_interpreter_curve.c
 *
 * Compile:
 *   gcc -O2 -Wall -Wextra -std=c99 test_interpreter_curve.c expr.c curve.c pdb.c -lm -o test_interpreter_curve
 *
 * Usage:
 *   ./test_interpreter_curve --t0 0 --t1 10 --n 50 --M 2000 --x "t" --y "0" --z "0" > line.pdb
 *   ./test_interpreter_curve --t0 0 --t1 18.8495559215 --n 200 --x "cos(t)" --y "sin(t)" --z "t" > helix.pdb
 *   ./test_interpreter_curve --t0 -2 --t1 2 --n 80 --x "t" --y "t*t" --z "0" > parabola.pdb
 *
 * Notes:
 * - Produces a PDB polyline: one CA atom per point.
 * - Uses curve_ws_build + curve_ws_resample_equal_arclen for equal arc-length output.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "expr.h"
#include "curve.h"

#ifndef ERRBUF_SZ
#define ERRBUF_SZ 256
#endif

typedef struct {
  const char *xs, *ys, *zs;
  Expr ex, ey, ez;
  char err[ERRBUF_SZ];
} ExprCurve3;

static int exprcurve3_init(ExprCurve3 *C, const char *xs, const char *ys, const char *zs){
  memset(C, 0, sizeof(*C));
  C->xs = xs; C->ys = ys; C->zs = zs;

  expr_init(&C->ex);
  expr_init(&C->ey);
  expr_init(&C->ez);

  ExprStatus st;

  st = expr_compile(&C->ex, xs);
  if (st != EXPR_OK){
    snprintf(C->err, sizeof(C->err), "compile x failed at %d: %s", C->ex.err_pos, expr_errmsg(&C->ex));
    return -1;
  }
  st = expr_compile(&C->ey, ys);
  if (st != EXPR_OK){
    snprintf(C->err, sizeof(C->err), "compile y failed at %d: %s", C->ey.err_pos, expr_errmsg(&C->ey));
    return -1;
  }
  st = expr_compile(&C->ez, zs);
  if (st != EXPR_OK){
    snprintf(C->err, sizeof(C->err), "compile z failed at %d: %s", C->ez.err_pos, expr_errmsg(&C->ez));
    return -1;
  }

  if (expr_uses_var_x(&C->ex) || expr_uses_var_y(&C->ex)){
    snprintf(C->err, sizeof(C->err), "x-expr must not reference x or y (only t allowed)");
    return -1;
  }
  if (expr_uses_var_y(&C->ey)){
    snprintf(C->err, sizeof(C->err), "y-expr must not reference y (self reference not allowed); use t and/or x");
    return -1;
  }
  /* z-expr may reference x and y: no restriction here */

  return 0;
}

static void exprcurve3_free(ExprCurve3 *C){
  expr_free(&C->ex);
  expr_free(&C->ey);
  expr_free(&C->ez);
}

static void die_eval(const char *comp, double t, const char *expr_str, Expr *E){
  fprintf(stderr, "eval error: component=%s t=%.17g expr=\"%s\": %s\n",
          comp, t, expr_str ? expr_str : "(null)", expr_errmsg(E));
  exit(1);
}

static Vec3 curve_pos_expr3(double t, void *user){
  ExprCurve3 *C = (ExprCurve3*)user;

  ExprEnv env;
  memset(&env, 0, sizeof(env));
  env.t = t;

  /* Evaluate in fixed order: x then y then z */
  if (expr_eval_env(&C->ex, &env, &env.x) != EXPR_OK || !isfinite(env.x))
    die_eval("x", t, C->xs, &C->ex);

  if (expr_eval_env(&C->ey, &env, &env.y) != EXPR_OK || !isfinite(env.y))
    die_eval("y", t, C->ys, &C->ey);

  double z = 0.0;
  if (expr_eval_env(&C->ez, &env, &z) != EXPR_OK || !isfinite(z))
    die_eval("z", t, C->zs, &C->ez);

  return (Vec3){ env.x, env.y, z };
}

static void usage(const char *argv0){
  fprintf(stderr,
    "Usage:\n"
    "  %s --x XEXPR --y YEXPR --z ZEXPR [--t0 T0] [--t1 T1] [--n N] [--M M]\n"
    "\n"
    "Output (to stdout): lines of 'C x y z'\n"
    "Expressions may use variables: t, x, y\n"
    "Dependency rules:\n"
    "  x-expr: may use t only\n"
    "  y-expr: may use t and x\n"
    "  z-expr: may use t, x, and y\n"
    "\n"
    "  --t0/--t1     parameter range (default 0..10)\n"
    "  --n           number of output points (equal arc-length) (default 100)\n"
    "  --M           oversample points for arclen table (default max(2000, 50*n))\n",
    argv0);
}

int main(int argc, char **argv){
  const char *xs = NULL, *ys = NULL, *zs = NULL;
  double t0 = 0.0, t1 = 10.0;
  int n = 100;
  int M = 0;

  for (int i=1;i<argc;i++){
    if (!strcmp(argv[i], "--x") && i+1<argc){ xs = argv[++i]; continue; }
    if (!strcmp(argv[i], "--y") && i+1<argc){ ys = argv[++i]; continue; }
    if (!strcmp(argv[i], "--z") && i+1<argc){ zs = argv[++i]; continue; }
    if (!strcmp(argv[i], "--t0") && i+1<argc){ t0 = atof(argv[++i]); continue; }
    if (!strcmp(argv[i], "--t1") && i+1<argc){ t1 = atof(argv[++i]); continue; }
    if (!strcmp(argv[i], "--n")  && i+1<argc){ n  = atoi(argv[++i]); continue; }
    if (!strcmp(argv[i], "--M")  && i+1<argc){ M  = atoi(argv[++i]); continue; }
    if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")){ usage(argv[0]); return 0; }

    fprintf(stderr, "Unknown/incomplete option: %s\n", argv[i]);
    usage(argv[0]);
    return 2;
  }

  if (!xs || !ys || !zs){
    fprintf(stderr, "Need --x, --y, --z\n");
    usage(argv[0]);
    return 2;
  }
  if (n < 2){
    fprintf(stderr, "--n must be >= 2\n");
    return 2;
  }
  if (t1 == t0){
    fprintf(stderr, "--t1 must differ from --t0\n");
    return 2;
  }
  if (M <= 0){
    M = (n < 2000) ? 2000 : (50*n);
  }
  if (M < 2) M = 2;

  ExprCurve3 C;
  if (exprcurve3_init(&C, xs, ys, zs) != 0){
    fprintf(stderr, "%s\n", C.err);
    exprcurve3_free(&C);
    return 1;
  }

  Vec3 *ws_p = (Vec3*)malloc((size_t)M * sizeof(Vec3));
  double *ws_s = (double*)malloc((size_t)M * sizeof(double));
  if (!ws_p || !ws_s){
    perror("malloc");
    free(ws_p); free(ws_s);
    exprcurve3_free(&C);
    return 1;
  }
  CurveWS ws = { .p = ws_p, .s = ws_s, .M = M };

  if (curve_ws_build(curve_pos_expr3, &C, t0, t1, &ws) != 0){
    fprintf(stderr, "curve_ws_build failed\n");
    free(ws_p); free(ws_s);
    exprcurve3_free(&C);
    return 1;
  }

  Vec3 *pts = (Vec3*)malloc((size_t)n * sizeof(Vec3));
  if (!pts){
    perror("malloc");
    free(ws_p); free(ws_s);
    exprcurve3_free(&C);
    return 1;
  }

  if (curve_ws_resample_equal_arclen(&ws, n, pts) != 0){
    fprintf(stderr, "curve_ws_resample_equal_arclen failed\n");
    free(pts);
    free(ws_p); free(ws_s);
    exprcurve3_free(&C);
    return 1;
  }

  for (int i=0;i<n;i++){
    printf("C %.17g %.17g %.17g\n", pts[i].x, pts[i].y, pts[i].z);
  }

  free(pts);
  free(ws_p); free(ws_s);
  exprcurve3_free(&C);
  return 0;
}
