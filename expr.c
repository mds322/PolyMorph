/* expr.c - minimal expression compiler/evaluator (infix -> RPN via shunting-yard)
 *Generated with significant support from GPT-5.2
 * Variables supported in expressions:
 *   t, x, y
 * Intended evaluation order for curve components:
 *   x = XEXPR(t)
 *   y = YEXPR(t, x)
 *   z = ZEXPR(t, x, y)
 * Notes:
 * - Compile produces E->rpn[] tokens (RPN program).
 * - Evaluation runs RPN with a fixed-size value stack.
 */
#include "expr.h"

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>

/* record a message + best-effort position */
static void seterr(Expr *E, const char *msg, int pos){
  if (!E) return;
  snprintf(E->err_msg, sizeof(E->err_msg), "%s", msg ? msg : "error");
  E->err_pos = pos;
}

const char *expr_errmsg(const Expr *E){
  return (E && E->err_msg[0]) ? E->err_msg : "error";
}

static char *xstrdup(const char *s){
  size_t n = strlen(s);
  char *p = (char*)malloc(n+1);
  if (!p) return NULL;
  memcpy(p, s, n+1);
  return p;
}

/* ---- function table ---- */

static int func_find(const Expr *E, const char *name){
  for (int i=0;i<E->nfunc;i++){
    if (strcmp(E->funcs[i].name, name) == 0) return i;
  }
  return -1;
}

static int func_add(Expr *E, const char *name, int arity, ExprFn1 f1, ExprFn2 f2){
  if (!E || !name || !*name) return -1;
  if (arity != 1 && arity != 2) return -1;
  if (func_find(E, name) >= 0) return -1;

  if (E->nfunc == E->capfunc){
    int newcap = (E->capfunc == 0) ? 32 : (E->capfunc * 2);
    ExprFunc *nf = (ExprFunc*)realloc(E->funcs, (size_t)newcap * sizeof(ExprFunc));
    if (!nf) return -1;
    E->funcs = nf;
    E->capfunc = newcap;
  }

  char *nm = xstrdup(name);
  if (!nm) return -1;

  ExprFunc *F = &E->funcs[E->nfunc++];
  F->name = nm;
  F->arity = arity;
  if (arity == 1) F->fn.f1 = f1;
  else            F->fn.f2 = f2;

  return 0;
}

static double fn_abs(double x){ return fabs(x); }
static double fn_min(double a,double b){ return (a<b)?a:b; }
static double fn_max(double a,double b){ return (a>b)?a:b; }

/* ---- public API ---- */

void expr_init(Expr *E){
  if (!E) return;
  memset(E, 0, sizeof(*E));
  E->err_msg[0] = '\0';
  E->err_pos = 0;

  /* trig */
  expr_register_fn1(E, "sin", sin);
  expr_register_fn1(E, "cos", cos);
  expr_register_fn1(E, "tan", tan);
  expr_register_fn1(E, "asin", asin);
  expr_register_fn1(E, "acos", acos);
  expr_register_fn1(E, "atan", atan);

  /* hyperbolic */
  expr_register_fn1(E, "sinh", sinh);
  expr_register_fn1(E, "cosh", cosh);
  expr_register_fn1(E, "tanh", tanh);
  expr_register_fn1(E, "asinh", asinh);
  expr_register_fn1(E, "acosh", acosh);
  expr_register_fn1(E, "atanh", atanh);

  /* misc */
  expr_register_fn1(E, "sqrt", sqrt);
  expr_register_fn1(E, "exp", exp);
  expr_register_fn1(E, "log", log);
  expr_register_fn1(E, "log10", log10);
  expr_register_fn1(E, "abs", fn_abs);
  expr_register_fn1(E, "floor", floor);
  expr_register_fn1(E, "ceil", ceil);

  /* 2-arg */
  expr_register_fn2(E, "pow", pow);
  expr_register_fn2(E, "atan2", atan2);
  expr_register_fn2(E, "min", fn_min);
  expr_register_fn2(E, "max", fn_max);
}

void expr_free(Expr *E){
  if (!E) return;

  free(E->rpn);
  E->rpn = NULL;
  E->nrpn = 0;

  for (int i=0;i<E->nfunc;i++){
    free((void*)E->funcs[i].name);
  }
  free(E->funcs);
  E->funcs = NULL;
  E->nfunc = 0;
  E->capfunc = 0;

  E->err_msg[0] = '\0';
  E->err_pos = 0;
}

int expr_register_fn1(Expr *E, const char *name, ExprFn1 f){
  return func_add(E, name, 1, f, NULL);
}

int expr_register_fn2(Expr *E, const char *name, ExprFn2 f){
  return func_add(E, name, 2, NULL, f);
}

/* ---- tokenizer + shunting-yard ---- */

typedef enum { TT_NUM, TT_IDENT, TT_OP, TT_LP, TT_RP, TT_COMMA, TT_END, TT_BAD } TokType;

typedef struct {
  TokType ty;
  double num;
  char ident[64];
  int op;    /* operator char, or offending char for TT_BAD */
  int pos;   /* byte offset in source string */
} Tok;

static int isident0(int c){ return isalpha((unsigned char)c) || c=='_'; }
static int isident(int c){ return isalnum((unsigned char)c) || c=='_'; }

static Tok next_tok(const char *expr_str, int *ioff){
  int i = *ioff;
  while (expr_str[i] && isspace((unsigned char)expr_str[i])) i++;

  Tok tok;
  memset(&tok, 0, sizeof(tok));
  tok.pos = i;

  if (!expr_str[i]){
    tok.ty = TT_END;
    *ioff = i;
    return tok;
  }

  char c = expr_str[i];

  if ((c>='0' && c<='9') || c=='.'){
    char *end = NULL;
    tok.num = strtod(expr_str+i, &end);
    tok.ty = TT_NUM;
    i = (int)(end - expr_str);
    *ioff = i;
    return tok;
  }

  if (isident0((unsigned char)c)){
    int j = 0;
    while (expr_str[i] && isident((unsigned char)expr_str[i])){
      if (j + 1 >= (int)sizeof(tok.ident)){
        while (expr_str[i] && isident((unsigned char)expr_str[i])) i++;
        tok.ident[0] = '\0';
        tok.ty = TT_IDENT;
        *ioff = i;
        return tok;
      }
      tok.ident[j++] = expr_str[i++];
    }
    tok.ident[j] = '\0';
    tok.ty = TT_IDENT;
    *ioff = i;
    return tok;
  }

  if (c=='+'||c=='-'||c=='*'||c=='/'||c=='^'){
    tok.ty = TT_OP;
    tok.op = c;
    *ioff = i+1;
    return tok;
  }
  if (c=='('){ tok.ty = TT_LP; *ioff = i+1; return tok; }
  if (c==')'){ tok.ty = TT_RP; *ioff = i+1; return tok; }
  if (c==','){ tok.ty = TT_COMMA; *ioff = i+1; return tok; }

  tok.ty = TT_BAD;
  tok.op = (unsigned char)c;
  *ioff = i+1;
  return tok;
}

static int op_prec(int op){
  switch (op){
    case '~': return 4; /* unary minus */
    case '^': return 3;
    case '*': case '/': return 2;
    case '+': case '-': return 1;
    default: return 0;
  }
}

static int op_right_assoc(int op){
  return (op == '^' || op == '~');
}

static int rpn_push(Expr *E, ExprTok tk, int *rpn_cap){
  if (E->nrpn == *rpn_cap){
    int newcap = (*rpn_cap==0) ? 64 : (*rpn_cap*2);
    ExprTok *nr = (ExprTok*)realloc(E->rpn, (size_t)newcap * sizeof(ExprTok));
    if (!nr) return -1;
    E->rpn = nr;
    *rpn_cap = newcap;
  }
  E->rpn[E->nrpn++] = tk;
  return 0;
}

static int opstack_push(int **op_stack, int *op_top, int *op_cap, int v){
  if (*op_top == *op_cap){
    int newcap = (*op_cap == 0) ? 64 : (*op_cap * 2);
    int *ns = (int*)realloc(*op_stack, (size_t)newcap * sizeof(int));
    if (!ns) return -1;
    *op_stack = ns;
    *op_cap = newcap;
  }
  (*op_stack)[(*op_top)++] = v;
  return 0;
}

ExprStatus expr_compile(Expr *E, const char *expr_str){
  if (!E || !expr_str) return EXPR_ERR_SYNTAX;

  free(E->rpn);
  E->rpn = NULL;
  E->nrpn = 0;
  E->err_msg[0] = '\0';
  E->err_pos = 0;

  int *op_stack = NULL;
  int op_cap = 0, op_top = 0;

  if (opstack_push(&op_stack, &op_top, &op_cap, 0) != 0){
    seterr(E, "alloc failed", 0);
    return EXPR_ERR_ALLOC;
  }
  op_top = 0;

  int rpn_cap = 0;
  int scan_pos = 0;

  Tok prev_tok;
  memset(&prev_tok, 0, sizeof(prev_tok));
  prev_tok.ty = TT_END;

  for (Tok tok = next_tok(expr_str, &scan_pos);
       tok.ty != TT_END;
       tok = next_tok(expr_str, &scan_pos)){

    if (tok.ty == TT_BAD){
      char msg[96];
      snprintf(msg, sizeof(msg), "unexpected character '%c'", (char)tok.op);
      free(op_stack);
      seterr(E, msg, tok.pos);
      return EXPR_ERR_SYNTAX;
    }

    if (tok.ty == TT_NUM){
      ExprTok insn = { .kind = EXRPN_NUM, .num = tok.num, .op = 0, .func_index = 0 };
      if (rpn_push(E, insn, &rpn_cap) != 0){
        free(op_stack);
        seterr(E, "alloc failed", tok.pos);
        return EXPR_ERR_ALLOC;
      }
      prev_tok = tok;
      continue;
    }

    if (tok.ty == TT_IDENT){
      if (tok.ident[0] == '\0'){
        free(op_stack);
        seterr(E, "identifier too long", tok.pos);
        return EXPR_ERR_SYNTAX;
      }

      /* variables */
      if (strcmp(tok.ident, "t") == 0){
        ExprTok insn = { .kind = EXRPN_VAR_T };
        if (rpn_push(E, insn, &rpn_cap) != 0){
          free(op_stack);
          seterr(E, "alloc failed", tok.pos);
          return EXPR_ERR_ALLOC;
        }
        prev_tok = tok;
        continue;
      }
      if (strcmp(tok.ident, "x") == 0){
        ExprTok insn = { .kind = EXRPN_VAR_X };
        if (rpn_push(E, insn, &rpn_cap) != 0){
          free(op_stack);
          seterr(E, "alloc failed", tok.pos);
          return EXPR_ERR_ALLOC;
        }
        prev_tok = tok;
        continue;
      }
      if (strcmp(tok.ident, "y") == 0){
        ExprTok insn = { .kind = EXRPN_VAR_Y };
        if (rpn_push(E, insn, &rpn_cap) != 0){
          free(op_stack);
          seterr(E, "alloc failed", tok.pos);
          return EXPR_ERR_ALLOC;
        }
        prev_tok = tok;
        continue;
      }

      /* constants */
      if (strcmp(tok.ident, "pi") == 0 || strcmp(tok.ident, "PI") == 0){
        ExprTok insn = { .kind = EXRPN_NUM, .num = M_PI };
        if (rpn_push(E, insn, &rpn_cap) != 0){
          free(op_stack);
          seterr(E, "alloc failed", tok.pos);
          return EXPR_ERR_ALLOC;
        }
        prev_tok = tok;
        continue;
      }
      if (strcmp(tok.ident, "e") == 0 || strcmp(tok.ident, "E") == 0){
        ExprTok insn = { .kind = EXRPN_NUM, .num = M_E };
        if (rpn_push(E, insn, &rpn_cap) != 0){
          free(op_stack);
          seterr(E, "alloc failed", tok.pos);
          return EXPR_ERR_ALLOC;
        }
        prev_tok = tok;
        continue;
      }

      /* function */
      int func_index = func_find(E, tok.ident);
      if (func_index < 0){
        free(op_stack);
        seterr(E, "unknown identifier", tok.pos);
        return EXPR_ERR_UNKNOWN_IDENT;
      }

      if (opstack_push(&op_stack, &op_top, &op_cap, -(func_index+1)) != 0){
        free(op_stack);
        seterr(E, "alloc failed", tok.pos);
        return EXPR_ERR_ALLOC;
      }
      prev_tok = tok;
      continue;
    }

    if (tok.ty == TT_OP){
      int op = tok.op;

      int prev_is_value = (prev_tok.ty == TT_NUM || prev_tok.ty == TT_IDENT || prev_tok.ty == TT_RP);
      if (!prev_is_value){
        if (op == '-') op = '~';
        else if (op == '+'){ prev_tok = tok; continue; }
      }

      while (op_top > 0){
        int stack_top = op_stack[op_top-1];
        if (stack_top == '(' || stack_top < 0) break;

        int p_top = op_prec(stack_top);
        int p_op  = op_prec(op);
        if (p_top > p_op || (p_top == p_op && !op_right_assoc(op))){
          op_top--;
          ExprTok insn = { .kind = EXRPN_OP, .op = stack_top };
          if (rpn_push(E, insn, &rpn_cap) != 0){
            free(op_stack);
            seterr(E, "alloc failed", tok.pos);
            return EXPR_ERR_ALLOC;
          }
        } else break;
      }

      if (opstack_push(&op_stack, &op_top, &op_cap, op) != 0){
        free(op_stack);
        seterr(E, "alloc failed", tok.pos);
        return EXPR_ERR_ALLOC;
      }

      prev_tok = tok;
      continue;
    }

    if (tok.ty == TT_LP){
      if (opstack_push(&op_stack, &op_top, &op_cap, '(') != 0){
        free(op_stack);
        seterr(E, "alloc failed", tok.pos);
        return EXPR_ERR_ALLOC;
      }
      prev_tok = tok;
      continue;
    }

    if (tok.ty == TT_COMMA){
      while (op_top > 0 && op_stack[op_top-1] != '('){
        int stack_top = op_stack[--op_top];
        if (stack_top < 0){
          free(op_stack);
          seterr(E, "misplaced comma", tok.pos);
          return EXPR_ERR_SYNTAX;
        }
        ExprTok insn = { .kind = EXRPN_OP, .op = stack_top };
        if (rpn_push(E, insn, &rpn_cap) != 0){
          free(op_stack);
          seterr(E, "alloc failed", tok.pos);
          return EXPR_ERR_ALLOC;
        }
      }
      if (op_top == 0){
        free(op_stack);
        seterr(E, "misplaced comma", tok.pos);
        return EXPR_ERR_SYNTAX;
      }
      prev_tok = tok;
      continue;
    }

    if (tok.ty == TT_RP){
      while (op_top > 0 && op_stack[op_top-1] != '('){
        int stack_top = op_stack[--op_top];
        if (stack_top < 0){
          free(op_stack);
          seterr(E, "bad ')'", tok.pos);
          return EXPR_ERR_SYNTAX;
        }
        ExprTok insn = { .kind = EXRPN_OP, .op = stack_top };
        if (rpn_push(E, insn, &rpn_cap) != 0){
          free(op_stack);
          seterr(E, "alloc failed", tok.pos);
          return EXPR_ERR_ALLOC;
        }
      }
      if (op_top == 0){
        free(op_stack);
        seterr(E, "unmatched ')'", tok.pos);
        return EXPR_ERR_SYNTAX;
      }
      op_top--; /* pop '(' */

      if (op_top > 0 && op_stack[op_top-1] < 0){
        int func_index = -op_stack[--op_top] - 1;
        ExprFunc *F = &E->funcs[func_index];

        ExprTok insn;
        memset(&insn, 0, sizeof(insn));
        insn.kind = (F->arity == 1) ? EXRPN_FN1 : EXRPN_FN2;
        insn.func_index = func_index;

        if (rpn_push(E, insn, &rpn_cap) != 0){
          free(op_stack);
          seterr(E, "alloc failed", tok.pos);
          return EXPR_ERR_ALLOC;
        }
      }

      prev_tok = tok;
      continue;
    }

    free(op_stack);
    seterr(E, "unexpected token", tok.pos);
    return EXPR_ERR_SYNTAX;
  }

  while (op_top > 0){
    int stack_top = op_stack[--op_top];
    if (stack_top == '('){
      free(op_stack);
      seterr(E, "unmatched '('", scan_pos);
      return EXPR_ERR_SYNTAX;
    }
    if (stack_top < 0){
      free(op_stack);
      seterr(E, "function missing ()", scan_pos);
      return EXPR_ERR_SYNTAX;
    }

    ExprTok insn = { .kind = EXRPN_OP, .op = stack_top };
    if (rpn_push(E, insn, &rpn_cap) != 0){
      free(op_stack);
      seterr(E, "alloc failed", scan_pos);
      return EXPR_ERR_ALLOC;
    }
  }

  free(op_stack);
  return EXPR_OK;
}

/* ---- helpers ---- */

int expr_uses_var_x(const Expr *E){
  if (!E) return 0;
  for (int i=0;i<E->nrpn;i++) if (E->rpn[i].kind == EXRPN_VAR_X) return 1;
  return 0;
}

int expr_uses_var_y(const Expr *E){
  if (!E) return 0;
  for (int i=0;i<E->nrpn;i++) if (E->rpn[i].kind == EXRPN_VAR_Y) return 1;
  return 0;
}

/* ---- evaluation ---- */

static ExprStatus eval_seterr(Expr *E, ExprStatus st, const char *msg){
  seterr(E, msg, 0);
  return st;
}

ExprStatus expr_eval_env(Expr *E, const ExprEnv *env, double *out){
  if (!E || !env || !out) return EXPR_ERR_SYNTAX;

  double val_stack[256];
  int val_top = 0;

  for (int i=0;i<E->nrpn;i++){
    ExprTok insn = E->rpn[i];

    if (insn.kind == EXRPN_NUM){
      if (val_top >= (int)(sizeof(val_stack)/sizeof(val_stack[0])))
        return eval_seterr(E, EXPR_ERR_STACK, "stack overflow");
      val_stack[val_top++] = insn.num;
      continue;
    }

    if (insn.kind == EXRPN_VAR_T){
      if (val_top >= (int)(sizeof(val_stack)/sizeof(val_stack[0])))
        return eval_seterr(E, EXPR_ERR_STACK, "stack overflow");
      val_stack[val_top++] = env->t;
      continue;
    }

    if (insn.kind == EXRPN_VAR_X){
      if (val_top >= (int)(sizeof(val_stack)/sizeof(val_stack[0])))
        return eval_seterr(E, EXPR_ERR_STACK, "stack overflow");
      val_stack[val_top++] = env->x;
      continue;
    }

    if (insn.kind == EXRPN_VAR_Y){
      if (val_top >= (int)(sizeof(val_stack)/sizeof(val_stack[0])))
        return eval_seterr(E, EXPR_ERR_STACK, "stack overflow");
      val_stack[val_top++] = env->y;
      continue;
    }

    if (insn.kind == EXRPN_OP){
      if (insn.op == '~'){
        if (val_top < 1) return eval_seterr(E, EXPR_ERR_STACK, "stack underflow");
        val_stack[val_top-1] = -val_stack[val_top-1];
        if (!isfinite(val_stack[val_top-1]))
          return eval_seterr(E, EXPR_ERR_DOMAIN, "non-finite (NaN/Inf)");
        continue;
      }

      if (val_top < 2) return eval_seterr(E, EXPR_ERR_STACK, "stack underflow");
      double b = val_stack[--val_top];
      double a = val_stack[--val_top];
      double r = 0.0;

      switch (insn.op){
        case '+': r = a + b; break;
        case '-': r = a - b; break;
        case '*': r = a * b; break;
        case '/': r = a / b; break;
        case '^': r = pow(a, b); break;
        default:  return eval_seterr(E, EXPR_ERR_SYNTAX, "bad operator");
      }

      if (!isfinite(r))
        return eval_seterr(E, EXPR_ERR_DOMAIN, "non-finite (NaN/Inf)");

      val_stack[val_top++] = r;
      continue;
    }

    if (insn.kind == EXRPN_FN1){
      if (val_top < 1) return eval_seterr(E, EXPR_ERR_STACK, "stack underflow");
      if (insn.func_index < 0 || insn.func_index >= E->nfunc)
        return eval_seterr(E, EXPR_ERR_SYNTAX, "bad function index");

      double a = val_stack[val_top-1];
      const ExprFunc *F = &E->funcs[insn.func_index];
      double r = F->fn.f1(a);

      if (!isfinite(r))
        return eval_seterr(E, EXPR_ERR_DOMAIN, "non-finite (NaN/Inf)");

      val_stack[val_top-1] = r;
      continue;
    }

    if (insn.kind == EXRPN_FN2){
      if (val_top < 2) return eval_seterr(E, EXPR_ERR_STACK, "stack underflow");
      if (insn.func_index < 0 || insn.func_index >= E->nfunc)
        return eval_seterr(E, EXPR_ERR_SYNTAX, "bad function index");

      double b = val_stack[--val_top];
      double a = val_stack[--val_top];
      const ExprFunc *F = &E->funcs[insn.func_index];
      double r = F->fn.f2(a, b);

      if (!isfinite(r))
        return eval_seterr(E, EXPR_ERR_DOMAIN, "non-finite (NaN/Inf)");

      val_stack[val_top++] = r;
      continue;
    }

    return eval_seterr(E, EXPR_ERR_SYNTAX, "bad token");
  }

  if (val_top != 1) return eval_seterr(E, EXPR_ERR_STACK, "stack not singular at end");

  *out = val_stack[0];
  return EXPR_OK;
}
