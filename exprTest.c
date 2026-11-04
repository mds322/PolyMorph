#include "expr.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

/* optional example "special function" */
static double fn_square(double x){ return x*x; }

static int expect_ok_eval(const char *s, double t, double want, double tol){
  Expr E;
  expr_init(&E);
  ExprStatus rc = expr_compile(&E, s);
  if (rc != EXPR_OK){
    fprintf(stderr, "FAIL compile: '%s' rc=%d msg=%s pos=%d\n", s, rc, expr_errmsg(&E), E.err_pos);
    expr_free(&E);
    return 1;
  }
  double out=0.0;
  rc = expr_eval(&E, t, &out);
  if (rc != EXPR_OK){
    fprintf(stderr, "FAIL eval: '%s' rc=%d\n", s, rc);
    expr_free(&E);
    return 1;
  }
  if (fabs(out - want) > tol){
    fprintf(stderr, "FAIL value: '%s' t=%g got=%g want=%g\n", s, t, out, want);
    expr_free(&E);
    return 1;
  }
  expr_free(&E);
  return 0;
}

static int expect_compile_fail(const char *s){
  Expr E;
  expr_init(&E);
  ExprStatus rc = expr_compile(&E, s);
  if (rc == EXPR_OK){
    fprintf(stderr, "FAIL expected compile error but compiled: '%s'\n", s);
    expr_free(&E);
    return 1;
  }
  expr_free(&E);
  return 0;
}

static int expect_eval_fail_stack(const char *s){
  Expr E;
  expr_init(&E);
  ExprStatus rc = expr_compile(&E, s);
  if (rc != EXPR_OK){
    fprintf(stderr, "FAIL compile (expected eval fail): '%s' rc=%d msg=%s\n", s, rc, expr_errmsg(&E));
    expr_free(&E);
    return 1;
  }
  double out=0.0;
  rc = expr_eval(&E, 1.0, &out);
  if (rc != EXPR_ERR_STACK){
    fprintf(stderr, "FAIL expected EXPR_ERR_STACK, got rc=%d for '%s'\n", rc, s);
    expr_free(&E);
    return 1;
  }
  expr_free(&E);
  return 0;
}

static int test_long_identifier_no_overflow(void){
  /* ident buffer in lexer is 64; should truncate, then become unknown identifier => compile fail */
  char s[512];
  memset(s, 'a', sizeof(s));
  s[0] = 'n'; s[1] = 'a'; s[2] = 'm'; s[3] = 'e';
  s[4] = '_';
  s[sizeof(s)-2] = '1';
  s[sizeof(s)-1] = '\0';
  return expect_compile_fail(s);
}

static int test_deep_stack(void){
  /* Build expression that pushes too many numbers before any ops.
     Our eval stack is 256. This should return EXPR_ERR_STACK, not crash. */
  char s[4096];
  s[0]='\0';
  for (int i=0;i<400;i++){
    strcat(s, "1 ");
    if (strlen(s) > sizeof(s)-10) break;
  }
  return expect_eval_fail_stack(s);
}

static int test_unknown_identifier(void){
  return expect_compile_fail("system(1)"); /* unknown identifier; should not execute */
}

static int test_register_function(void){
  Expr E;
  expr_init(&E);
  if (expr_register_fn1(&E, "square", fn_square) != 0){
    fprintf(stderr, "FAIL registering square()\n");
    expr_free(&E);
    return 1;
  }
  if (expr_compile(&E, "square(t) + square(2)") != EXPR_OK){
    fprintf(stderr, "FAIL compile square expr: %s\n", expr_errmsg(&E));
    expr_free(&E);
    return 1;
  }
  double out=0.0;
  if (expr_eval(&E, 3.0, &out) != EXPR_OK){
    fprintf(stderr, "FAIL eval square expr\n");
    expr_free(&E);
    return 1;
  }
  if (fabs(out - (9.0 + 4.0)) > 1e-12){
    fprintf(stderr, "FAIL square value got=%g\n", out);
    expr_free(&E);
    return 1;
  }
  expr_free(&E);
  return 0;
}

int main(void){
  int fails = 0;

  fails += expect_ok_eval("1+2*3", 0.0, 7.0, 1e-12);
  fails += expect_ok_eval("sin(pi/2)", 0.0, 1.0, 1e-12);
  fails += expect_ok_eval("cosh(0)", 0.0, 1.0, 1e-12);
  fails += expect_ok_eval("atan2(1,1)", 0.0, M_PI/4.0, 1e-12);
  fails += expect_ok_eval("-t*t + 2*t", 3.0, -9.0 + 6.0, 1e-12);
  fails += expect_ok_eval("pow(2,3)+2^3", 0.0, 16.0, 1e-12);

  fails += expect_compile_fail("sin)");          /* syntax */
  fails += expect_compile_fail("unknown(t)");    /* unknown identifier */
  fails += expect_compile_fail("sin(t");         /* unmatched '(' */
  fails += test_long_identifier_no_overflow();
  fails += test_unknown_identifier();
  fails += test_deep_stack();
  fails += test_register_function();

  if (fails == 0){
    printf("expr_sanity: OK\n");
    return 0;
  } else {
    printf("expr_sanity: FAIL (%d)\n", fails);
    return 1;
  }
}
