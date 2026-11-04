#ifndef CLI_H
#define CLI_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  const char *infile;
  const char *outfile;

  int k;
  int n_points;
  double t0, t1;

  const char *sel1_str;
  const char *sel2_str;

  /* curve expr strings */
  const char *xexpr;
  const char *yexpr;
  const char *zexpr;

/* Optional auto-tuning flag value from argv --t1 auto */
  int t1_auto;

/* Bond length options (possible replacement for auto scaling fixes*/
  double bond_len;

} CliOpt;


/* function that prints usage */

void cli_usage(const char *argv0);

/* Parse argv into CliOpt.
   Returns 0 on success, nonzero on error (prints message to stderr). */
int cli_parse(int argc, char **argv, CliOpt *opt);

#ifdef __cplusplus
}
#endif
#endif
