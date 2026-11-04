#ifndef CURVE_EXPR3_H
#define CURVE_EXPR3_H

#include "curve.h"
#include "expr.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  const char *xs, *ys, *zs;
  Expr ex, ey, ez;
} ExprCurve3;

static inline void exprcurve3_destroy(ExprCurve3 *C){
  if (!C) return;
  expr_free(&C->ex);
  expr_free(&C->ey);
  expr_free(&C->ez);
  memset(C, 0, sizeof(*C));
}

static inline void exprcurve3_die_eval(const char *comp, double t,
                                       const char *expr_str, Expr *E)
{
  fprintf(stderr, "curve expr error: component=%s t=%.17g expr=\"%s\": %s\n",
          comp, t, expr_str ? expr_str : "(null)", expr_errmsg(E));
  exit(1);
}

/* Compile/validate three expressions (Option-1 dependency rules).
   Returns 0 on success, nonzero on failure (prints message to stderr). */
static inline int exprcurve3_init(ExprCurve3 *C,
                                  const char *xs, const char *ys, const char *zs)
{
  if (!C || !xs || !ys || !zs) return -1;

  memset(C, 0, sizeof(*C));
  C->xs = xs; C->ys = ys; C->zs = zs;

  expr_init(&C->ex);
  expr_init(&C->ey);
  expr_init(&C->ez);

  if (expr_compile(&C->ex, xs) != EXPR_OK){
    fprintf(stderr, "compile x failed at %d: %s\n", C->ex.err_pos, expr_errmsg(&C->ex));
    exprcurve3_destroy(C);
    return -1;
  }
  if (expr_compile(&C->ey, ys) != EXPR_OK){
    fprintf(stderr, "compile y failed at %d: %s\n", C->ey.err_pos, expr_errmsg(&C->ey));
    exprcurve3_destroy(C);
    return -1;
  }
  if (expr_compile(&C->ez, zs) != EXPR_OK){
    fprintf(stderr, "compile z failed at %d: %s\n", C->ez.err_pos, expr_errmsg(&C->ez));
    exprcurve3_destroy(C);
    return -1;
  }

  /* Option-1 dependency rules */
  if (expr_uses_var_x(&C->ex) || expr_uses_var_y(&C->ex)){
    fprintf(stderr, "x-expr must not reference x or y (only t allowed)\n");
    exprcurve3_destroy(C);
    return -1;
  }
  if (expr_uses_var_y(&C->ey)){
    fprintf(stderr, "y-expr must not reference y (self reference not allowed); use t and/or x\n");
    exprcurve3_destroy(C);
    return -1;
  }

  return 0;
}

/* CurvePosFn adapter: evaluates x/y/z at parameter t.
   On eval error prints detailed error and exits(1) because CurvePosFn has no error return. */
static inline Vec3 curve_pos_expr3(double t, void *user){
  ExprCurve3 *C = (ExprCurve3*)user;

  ExprEnv env;
  memset(&env, 0, sizeof(env));
  env.t = t;

  if (expr_eval_env(&C->ex, &env, &env.x) != EXPR_OK || !isfinite(env.x))
    exprcurve3_die_eval("x", t, C->xs, &C->ex);

  if (expr_eval_env(&C->ey, &env, &env.y) != EXPR_OK || !isfinite(env.y))
    exprcurve3_die_eval("y", t, C->ys, &C->ey);

  double z = 0.0;
  if (expr_eval_env(&C->ez, &env, &z) != EXPR_OK || !isfinite(z))
    exprcurve3_die_eval("z", t, C->zs, &C->ez);

  return (Vec3){ env.x, env.y, z };
}

#ifdef __cplusplus
}
#endif
#endif
