#include "cli.h"
#include "vecops.h"
#include "curve.h"
#include "expr.h"
#include "curve_expr3.h"
#include "selection.h"
#include "scaling.h"
#include "pdb.h"
#include "curve_adapt.h"
#include "curvature.h"
#include "curve_scale.h"
#include "place.h"
#include "steric_scale.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

typedef struct {
  const Vec3 *unit_rel;
  size_t nAtoms;
  Vec3 axis_unit;
  PlaceOpt popt;
} PlaceUser;

static int place_cb(const Vec3 *pts, int K, void *user, Vec3 *out_coords){
  PlaceUser *U = (PlaceUser*)user;
  return place_parallel_transport(pts, K, U->unit_rel, U->nAtoms,
                                  U->axis_unit, &U->popt, out_coords);
}

static void dump_pts(const char *tag, const Vec3 *pts, int K){
  fprintf(stderr, "---- %s ----\n", tag);
  for (int i=0;i<K;i++){
    fprintf(stderr, "pts[%2d] = %.6f %.6f %.6f\n", i, pts[i].x, pts[i].y, pts[i].z);
  }
  for (int i=0;i<K-1;i++){
    fprintf(stderr, "step[%2d] = %.6f\n", i, vec3_dist(pts[i+1], pts[i]));
  }
}

static int dump_pts_pdb(const char *fname, const Vec3 *pts, int K){
  FILE *fp = fopen(fname, "w");
  if (!fp) return -1;
  for (int i=0;i<K;i++){
    fprintf(fp,
      "ATOM  %5d  PTS PTS A%4d    %8.3f%8.3f%8.3f  1.00  0.00\n",
      i+1, i+1, pts[i].x, pts[i].y, pts[i].z);
  }
  fprintf(fp, "END\n");
  fclose(fp);
  return 0;
}

/*Temporary location of bond-length buffering code...will move when testing complete*/
static double monomer_extent_along_axis(const Vec3 *unit_rel, size_t nAtoms, Vec3 axis_unit){
  if (!unit_rel || nAtoms == 0) return 0.0;
  axis_unit = vec3_nrm(axis_unit);
  double minp = 0.0, maxp = 0.0;
  for (size_t i=0;i<nAtoms;i++){
    double p = vec3_dot(unit_rel[i], axis_unit);
    if (i == 0){ minp = maxp = p; }
    else {
      if (p < minp) minp = p;
      if (p > maxp) maxp = p;
    }
  }
  return maxp - minp;
}


int main(int argc, char **argv){
  CliOpt opt;
  if (cli_parse(argc, argv, &opt) != 0) return 1;

  SelectionSpec s1, s2;
  if (selection_parse(opt.sel1_str, &s1) != 0){
    fprintf(stderr, "Bad --sel1: %s\n", opt.sel1_str);
    return 1;
  }
  if (selection_parse(opt.sel2_str, &s2) != 0){
    fprintf(stderr, "Bad --sel2: %s\n", opt.sel2_str);
    return 1;
  }

  int *line1 = (int*)malloc((size_t)PDB_MAX_ATOMS * sizeof(int));
  if (!line1){ perror("malloc line1"); return 1; }

  Atom *atoms = (Atom*)malloc((size_t)PDB_MAX_ATOMS * sizeof(Atom));
  if (!atoms){ perror("malloc atoms"); free(line1); return 1; }

  size_t nAtoms = 0;
  if (pdb_read_with_lines(opt.infile, atoms, PDB_MAX_ATOMS, &nAtoms, line1) != 0){
    free(atoms); free(line1);
    return 1;
  }
  if (nAtoms == 0){
    fprintf(stderr, "No ATOM/HETATM records read.\n");
    free(atoms); free(line1);
    return 1;
  }

  SelectionHit h1, h2;
  if (selection_resolve_pair(atoms, nAtoms, line1, &s1, &s2, &h1, &h2) != 0){
    fprintf(stderr, "Selection resolve failed for --sel1 \"%s\" --sel2 \"%s\".\n",
            opt.sel1_str, opt.sel2_str);
    free(atoms); free(line1);
    return 1;
  }

  Vec3 A = h1.atom.coord;
  Vec3 B = h2.atom.coord;
  Vec3 axis_unit = vec3_nrm(vec3_sub(B, A));
  Vec3 anchor    = vec3_scl(vec3_add(A, B), 0.5);

  int K = opt.k;

  Vec3 *unit_rel = (Vec3*)malloc(nAtoms * sizeof(Vec3));
  if (!unit_rel){
    perror("malloc unit_rel");
    free(atoms); free(line1);
    return 1;
  }
  for (size_t i=0;i<nAtoms;i++) unit_rel[i] = vec3_sub(atoms[i].coord, anchor);

  Vec3 *pts = (Vec3*)malloc((size_t)K * sizeof(Vec3));
  if (!pts){
    perror("malloc pts");
    free(atoms); free(unit_rel); free(line1);
    return 1;
  }


  double d_target = scaling_target_step_from_hits(&h1, &h2);
  if (d_target <= 0.0){
    fprintf(stderr, "Could not compute a valid distance from selection; defaulting to 5.0A.\n");
    d_target = 5.0;
  }

  double extent = monomer_extent_along_axis(unit_rel, nAtoms, axis_unit);
  double step_target = extent + opt.bond_len;

  /* fallback if something degenerate happened */
  if (!(step_target > 0.0)) step_target = d_target;

  fprintf(stderr,"spacing: d_sel=%.3fA extent_axis=%.3fA bond_len=%.3fA => step_target=%.3fA\n",
		  d_target, extent, opt.bond_len, step_target);

  int N_QUALITY = opt.n_points;
  if (N_QUALITY < 50) N_QUALITY = 50;

  int M = (N_QUALITY < 2000) ? 2000 : (N_QUALITY * 50);
  if (M < 2) M = 2;

  Vec3  *ws_p = (Vec3*)malloc((size_t)M * sizeof(Vec3));
  double*ws_s = (double*)malloc((size_t)M * sizeof(double));
  if (!ws_p || !ws_s){
    perror("malloc ws");
    free(ws_p); free(ws_s);
    free(atoms); free(pts); free(unit_rel); free(line1);
    return 1;
  }

  ExprCurve3 EC;
  if (exprcurve3_init(&EC, opt.xexpr, opt.yexpr, opt.zexpr) != 0){
    exprcurve3_destroy(&EC);
    free(ws_p); free(ws_s);
    free(atoms); free(pts); free(unit_rel); free(line1);
    return 1;
  }

  /* t1 auto / warnings */
  {
    CurveWS ws_tmp = { .p = ws_p, .s = ws_s, .M = M };
    const double s_goal = 0.25;

    if (opt.t1_auto){
      double target_len = (K > 1) ? ((double)(K - 1) * d_target) / s_goal : 0.0;
      double t1_new = opt.t1;
      if (t1_new == opt.t0) t1_new = opt.t0 + 1.0;

      if (curve_autotune_t1_for_target_len(curve_pos_expr3, &EC,
                                           opt.t0, &t1_new,
                                           target_len,
                                           &ws_tmp,
                                           32, 64, 1e-3) != 0){
        fprintf(stderr, "auto-t1 failed; try providing numeric --t1.\n");
        exprcurve3_destroy(&EC);
        free(ws_p); free(ws_s);
        free(atoms); free(pts); free(unit_rel); free(line1);
        return 1;
      }
      fprintf(stderr, "auto-t1: adjusted t1 from %.6g to %.6g (target scale %.3g)\n",
              opt.t1, t1_new, s_goal);
      opt.t1 = t1_new;
    } else {
      double L = curve_ws_arclen(curve_pos_expr3, &EC, opt.t0, opt.t1, &ws_tmp);
      if (L > 0.0 && K > 1){
        double s_needed = ((double)(K - 1) * d_target) / L;
        if (s_needed > 0.4){
          fprintf(stderr,
            "Warning: implied curve_scale=%.3g (>0.4). Consider --t1 auto or increasing |t1-t0|.\n",
            s_needed);
        }
      }
    }
  }

  /* adaptive build */
  CurveWS ws = { .p = ws_p, .s = ws_s, .M = M };
  enum { IT_MAX = 20 };
  const int MAX_M = 2000000;
  const double THETA_LIMIT_RAD = 10.0 * M_PI / 180.0;

  int ok = 0;
  for (int it=0; it<IT_MAX; it++){
    int rc = curve_ws_build_adaptive_turn(curve_pos_expr3, &EC, opt.t0, opt.t1,
                                          &ws, MAX_M, THETA_LIMIT_RAD);
    if (rc == 0){ ok = 1; break; }
    if (rc < 0){ ok = 0; break; }

    Vec3 *np = (Vec3*)realloc(ws.p, (size_t)ws.M * sizeof(Vec3));
    if (!np){ perror("realloc ws.p"); ok = 0; break; }
    ws.p = np;

    double *ns = (double*)realloc(ws.s, (size_t)ws.M * sizeof(double));
    if (!ns){ perror("realloc ws.s"); ok = 0; break; }
    ws.s = ns;
  }

  if (!ok){
    fprintf(stderr, "Adaptive curve build failed.\n");
    exprcurve3_destroy(&EC);
    free(ws.p); free(ws.s);
    free(atoms); free(pts); free(unit_rel); free(line1);
    return 1;
  }

  if (curve_ws_resample_equal_arclen(&ws, K, pts) != 0){
    fprintf(stderr, "Curve resample failed.\n");
    exprcurve3_destroy(&EC);
    free(ws.p); free(ws.s);
    free(atoms); free(pts); free(unit_rel); free(line1);
    return 1;
  }

  exprcurve3_destroy(&EC);
  free(ws.p);
  free(ws.s);

  /* DEBUG: inspect raw resampled points */
  dump_pts("pts_raw (after equal-arclen resample, before scaling)", pts, K);
  (void)dump_pts_pdb("pts_raw.pdb", pts, K);

  /* choose scale (step + curvature min) */
  CurveScaleReport rep;
  double scale0 = curve_choose_global_scale(pts, K, unit_rel, nAtoms, step_target, 3.0, &rep);
//  double scale0 = curve_choose_global_scale(pts, K, unit_rel, nAtoms, d_target, 3.0, &rep);
  if (!(scale0 > 0.0)) scale0 = 1.0;

  fprintf(stderr,
    "scale0 factors: step=%.6g curv=%.6g chosen=%.6g (Rmin=%.6g r_unit=%.6g)\n",
    rep.scale_step, rep.scale_curv, rep.scale_chosen, rep.Rmin, rep.r_unit);

  curve_scale_pts(pts, K, scale0);

  /* DEBUG: inspect scaled points */
  dump_pts("pts_scaled0 (after step+curv scaling)", pts, K);
  (void)dump_pts_pdb("pts_scaled0.pdb", pts, K);

  if ((size_t)K > 0 && nAtoms > SIZE_MAX / (size_t)K){
    fprintf(stderr, "Overflow: K*nAtoms too large\n");
    free(atoms); free(pts); free(unit_rel); free(line1);
    return 1;
  }
  size_t totalN = (size_t)K * nAtoms;

  Vec3 *out = (Vec3*)malloc(totalN * sizeof(Vec3));
  if (!out){
    perror("malloc out");
    free(atoms); free(pts); free(unit_rel); free(line1);
    return 1;
  }

  /* Prepare steric report even if steric loop is disabled */
  StericScaleReport srep;
  memset(&srep, 0, sizeof(srep));
  srep.final_scale_mult = 1.0;
  srep.last_dmin = HUGE_VAL;
  srep.converged = 0;

  PlaceUser pu;
  pu.unit_rel = unit_rel;
  pu.nAtoms = nAtoms;
  pu.axis_unit = axis_unit;
  pu.popt.twist_per_repeat = 0.0;

  /* ---- DEBUG STEP: place once (no steric loop) ---- */
  if (place_cb(pts, K, &pu, out) != 0){
    fprintf(stderr, "Placement failed.\n");
    free(atoms); free(out); free(pts); free(unit_rel); free(line1);
    return 1;
  }

  /* If/when enabling steric scaling, replace the one-shot place above with:
     StericScaleOpt sopt = { .W=6, .cutoff=1.2, .grow=1.03, .max_iter=20 };
     if (steric_scale_apply(pts, K, atoms, nAtoms, place_cb, &pu, out, &sopt, &srep) != 0) { ... }
     fprintf(stderr, "steric_scale: iters=%d mult=%.6g dmin=%.3f converged=%d\n", ...);
   */

  double scale_final = rep.scale_chosen * srep.final_scale_mult;

  Atom *meta_all = (Atom*)malloc(totalN * sizeof(Atom));
  if (!meta_all){
    perror("malloc meta_all");
    free(atoms); free(out); free(pts); free(unit_rel); free(line1);
    return 1;
  }

  int serial = 1;
  int resSeq0 = atoms[0].resSeq;
  int res_per_unit = pdb_count_residues_ordered(atoms, nAtoms);
  if (res_per_unit < 1) res_per_unit = 1;

  for (int r=0; r<K; r++){
    int res_off = r * res_per_unit;
    for (size_t i=0;i<nAtoms;i++){
      Atom a = atoms[i];
      a.serial = serial++;
      a.resSeq = resSeq0 + res_off + (a.resSeq - resSeq0);
      meta_all[(size_t)r*nAtoms + i] = a;
    }
  }

  if (pdb_write_atoms_with_conect_repeats(opt.outfile, meta_all, out, totalN,
                                         nAtoms, K, h1.index0, h2.index0) != 0){
    free(atoms); free(meta_all); free(out); free(pts); free(unit_rel); free(line1);
    return 1;
  }

  fprintf(stdout,
    "Wrote: repeats=%d unitAtoms=%zu totalAtoms=%zu quality=%d\n"
    "Scale: target_step=%.3fA curve_scale=%.6g (step=%.6g curv=%.6g steric_mult=%.6g) t0=%.6g t1=%.6g\n",
    K, nAtoms, totalN, N_QUALITY,
    d_target, scale_final, rep.scale_step, rep.scale_curv, srep.final_scale_mult,
    opt.t0, opt.t1);

  free(atoms);
  free(meta_all);
  free(out);
  free(pts);
  free(unit_rel);
  free(line1);
  return 0;
}
