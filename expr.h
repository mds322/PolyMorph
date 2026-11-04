#ifndef EXPR_H
#define EXPR_H

#include <stddef.h>

#ifndef M_PI
#define M_PI (3.14159265358979323846264338327950)
#endif /* M_PI */

#ifndef M_E
#define M_E  (2.71828182845904523536028747135266)
#endif /* M_E */

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
  EXPR_OK = 0,
  EXPR_ERR_SYNTAX = -1,
  EXPR_ERR_UNKNOWN_IDENT = -2,
  EXPR_ERR_STACK = -3,
  EXPR_ERR_ALLOC = -4,
  EXPR_ERR_DOMAIN = -5
} ExprStatus;

/* Function pointer types for runtime-registered functions */
typedef double (*ExprFn1)(double);
typedef double (*ExprFn2)(double, double);

typedef struct {
  const char *name;
  int arity; /* 1 or 2 */
  union { ExprFn1 f1; ExprFn2 f2; } fn;
} ExprFunc;

typedef enum {
  EXRPN_NUM   = 1,
  EXRPN_VAR_T = 2,
  EXRPN_VAR_X = 6,
  EXRPN_VAR_Y = 7,
  EXRPN_OP    = 3,
  EXRPN_FN1   = 4,
  EXRPN_FN2   = 5
} ExprRpnKind;

typedef struct {
  ExprRpnKind kind;
  double num;     /* for EXRPN_NUM */
  int op;         /* for EXRPN_OP: '+','-','*','/','^','~'(unary minus) */
  int func_index; /* for FN1/FN2: index into function table */
} ExprTok;

typedef struct {
  double t, x, y;
} ExprEnv;

typedef struct {
  ExprTok *rpn;
  int nrpn;

  /* function table (builtins + registered) */
  ExprFunc *funcs;
  int nfunc;
  int capfunc;

  char err_msg[128];
  int err_pos; /* byte offset in input string, best-effort */
} Expr;

/* Initialize expression engine with built-in functions. */
void expr_init(Expr *E);

/* Free memory owned by Expr. */
void expr_free(Expr *E);

/* Register a new function at runtime. Name is copied (so caller can free). */
int expr_register_fn1(Expr *E, const char *name, ExprFn1 f);
int expr_register_fn2(Expr *E, const char *name, ExprFn2 f);

/* Compile expression string into RPN. Overwrites previous compiled expression. */
ExprStatus expr_compile(Expr *E, const char *s);

/* Evaluate compiled expression for environment variables t/x/y. */
ExprStatus expr_eval_env(Expr *E, const ExprEnv *env, double *out);

/* Access last error message (also stored in E->err_msg). */
const char *expr_errmsg(const Expr *E);

/* Helper: check whether compiled expression references variable x or y */
int expr_uses_var_x(const Expr *E);
int expr_uses_var_y(const Expr *E);

#ifdef __cplusplus
}
#endif
#endif
