#include "cli.h"
#include "expr.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

static int streq(const char *a, const char *b){ return strcmp(a,b)==0; }

void cli_usage(const char *argv0){
  fprintf(stderr,
    "Usage:\n"
    "  %s --in IN.pdb --out OUT.pdb --k K --bond-len \\\n"
    "     --sel1 SPEC --sel2 SPEC \\\n"
    "     --curve expr XEXPR YEXPR ZEXPR\n"
    "\n"
    "Required:\n"
    "  --in FILE         input PDB\n"
    "  --out FILE        output PDB\n"
    "  --k K             number of repeats\n"
    "  --sel1 SPEC       selection for atom 1 (unit axis/linkage)\n"
    "  --sel2 SPEC       selection for atom 2 (unit axis/linkage)\n"
    "  --curve expr ...  parametric curve (expressions in quotes)\n"
    "\n"
    "Optional:\n"
    "  --n-points N      curve samples used for placement (default max(k,20))\n"
    "  --t0 T0           curve parameter start (default 0)\n"
    "  --t1 T1           curve parameter end   (default 10)\n"
    "\n"
    "Expr variables:\n"
    "  x(t)=XEXPR uses: t only\n"
    "  y(t)=YEXPR uses: t and x\n"
    "  z(t)=ZEXPR uses: t, x, y\n"
    , argv0);
}

static int need_arg(const char *opt, int i, int argc){
  if (i+1 >= argc){
    fprintf(stderr, "Error: %s requires an argument\n", opt);
    return 1;
  }
  return 0;
}

int cli_parse(int argc, char **argv, CliOpt *opt){
  memset(opt, 0, sizeof(*opt));
/*Set defaults*/
  opt->k = -1;
  opt->n_points = 0;
  opt->t0 = 0.0;
  opt->t1 = 10.0;
  opt->t1_auto =0;
  opt->bond_len = 1.8;

  for (int i=1;i<argc;i++){
    const char *a = argv[i];

    if (streq(a, "-h") || streq(a, "--help")){
      cli_usage(argv[0]);
      return 1;
    }

    if (streq(a, "--in")){
      if (need_arg(a,i,argc)) return 1;
      opt->infile = argv[++i];
      continue;
    }
    if (streq(a, "--out")){
      if (need_arg(a,i,argc)) return 1;
      opt->outfile = argv[++i];
      continue;
    }
    if (streq(a, "--k")){
      if (need_arg(a,i,argc)) return 1;
      opt->k = atoi(argv[++i]);
      continue;
    }
    if (streq(a, "--n-points") || streq(a, "--n_points")){
      if (need_arg(a,i,argc)) return 1;
      opt->n_points = atoi(argv[++i]);
      continue;
    }
    if (streq(a, "--t0")){
      if (need_arg(a,i,argc)) return 1;
      opt->t0 = atof(argv[++i]);
      continue;
    }
    if (streq(a, "--t1")){
      if (need_arg(a,i,argc)) return 1;
      const char *v = argv[++i];
      if(!strcmp(v, "auto")){
	      opt->t1_auto = 1;
      } else {
	      opt->t1 = atof(v);
	      opt->t1_auto = 0;
      }
      continue;
    }
    if (streq(a, "--sel1")){
      if (need_arg(a,i,argc)) return 1;
      opt->sel1_str = argv[++i];
      continue;
    }
    if (streq(a, "--sel2")){
      if (need_arg(a,i,argc)) return 1;
      opt->sel2_str = argv[++i];
      continue;
    }

    if (streq(a, "--curve")){
      if (need_arg(a,i,argc)) return 1;
      const char *kind = argv[++i];
      if (!streq(kind, "expr")){
        fprintf(stderr, "Error: only supported curve kind is: expr\n");
        return 1;
      }
      if (i + 3 >= argc){
        fprintf(stderr, "Error: --curve expr requires 3 expressions: XEXPR YEXPR ZEXPR\n");
        return 1;
      }
      opt->xexpr = argv[++i];
      opt->yexpr = argv[++i];
      opt->zexpr = argv[++i];
      continue;
    }

    if (streq(a, "--bond-len")){
      if (need_arg(a,i,argc)) return 1;
      opt->bond_len = atof(argv[++i]);
      continue;
    }

    fprintf(stderr, "Error: unknown option: %s\n", a);
    return 1;
  }

  /* required checks */
  if (!opt->infile){ fprintf(stderr, "Error: missing --in\n"); return 1; }
  if (!opt->outfile){ fprintf(stderr, "Error: missing --out\n"); return 1; }
  if (opt->k < 1){ fprintf(stderr, "Error: --k must be >= 1\n"); return 1; }
  if (!opt->sel1_str || !opt->sel2_str){
    fprintf(stderr, "Error: missing --sel1/--sel2\n");
    return 1;
  }
  if (!opt->xexpr || !opt->yexpr || !opt->zexpr){
    fprintf(stderr, "Error: missing --curve expr XEXPR YEXPR ZEXPR\n");
    return 1;
  }
  if (opt->t1 == opt->t0){
    fprintf(stderr, "Error: --t1 must differ from --t0\n");
    return 1;
  }
  if (opt->n_points < 2){
    opt->n_points = (opt->k > 20) ? opt->k : 20;
  }

  /*Check bond length was read from command line*/
  if (!(opt->bond_len >= 0.0)){
    fprintf(stderr, "Error: --bond-len must be >= 0\n");
    return 1;
  }

  /* Validate dependency rules by compiling quickly here.
     This ensures early user-friendly errors, before reading PDB etc. */
  {
    Expr ex, ey, ez;
    expr_init(&ex); expr_init(&ey); expr_init(&ez);

    if (expr_compile(&ex, opt->xexpr) != EXPR_OK){
      fprintf(stderr, "Error: XEXPR compile failed at %d: %s\n", ex.err_pos, expr_errmsg(&ex));
      expr_free(&ex); expr_free(&ey); expr_free(&ez);
      return 1;
    }
    if (expr_compile(&ey, opt->yexpr) != EXPR_OK){
      fprintf(stderr, "Error: YEXPR compile failed at %d: %s\n", ey.err_pos, expr_errmsg(&ey));
      expr_free(&ex); expr_free(&ey); expr_free(&ez);
      return 1;
    }
    if (expr_compile(&ez, opt->zexpr) != EXPR_OK){
      fprintf(stderr, "Error: ZEXPR compile failed at %d: %s\n", ez.err_pos, expr_errmsg(&ez));
      expr_free(&ex); expr_free(&ey); expr_free(&ez);
      return 1;
    }

    if (expr_uses_var_x(&ex) || expr_uses_var_y(&ex)){
      fprintf(stderr, "Error: XEXPR must not reference x or y (only t allowed)\n");
      expr_free(&ex); expr_free(&ey); expr_free(&ez);
      return 1;
    }
    if (expr_uses_var_y(&ey)){
      fprintf(stderr, "Error: YEXPR must not reference y (self reference not allowed); use t and/or x\n");
      expr_free(&ex); expr_free(&ey); expr_free(&ez);
      return 1;
    }

    expr_free(&ex); expr_free(&ey); expr_free(&ez);
  }

  return 0;
}
