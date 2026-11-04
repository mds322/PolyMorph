#include "steric_scale.h"
#include "clash_sortgrid.h"
#include "curve.h" /* curve_scale_pts */
#include <stdlib.h>
#include <string.h>
#include <math.h>

int steric_scale_apply(Vec3 *pts, int K,
                       const Atom *atoms, size_t nAtoms,
                       StericPlaceFn place_fn, void *place_user,
                       Vec3 *out_coords,
                       const StericScaleOpt *opt,
                       StericScaleReport *rep)
{
  if (!pts || K < 1) return -1;
  if (!atoms || nAtoms == 0) return -1;
  if (!place_fn || !out_coords) return -1;

  StericScaleOpt O;
  O.W = 6;
  O.cutoff = 1.2;
  O.grow = 1.03;
  O.max_iter = 20;
  if (opt) O = *opt;

  if (O.W < 1) O.W = 1;
  if (O.W > K-1) O.W = (K > 1) ? (K-1) : 1;
  if (!(O.cutoff > 0.0)) O.cutoff = 1.2;
  if (!(O.grow > 1.0)) O.grow = 1.03;
  if (O.max_iter < 0) O.max_iter = 0;

  StericScaleReport R;
  memset(&R, 0, sizeof(R));
  R.iters = 0;
  R.final_scale_mult = 1.0;
  R.last_dmin = HUGE_VAL;
  R.converged = 0;

  /* build keep index list once */
  int *idx_keep = (int*)malloc(nAtoms * sizeof(int));
  if (!idx_keep) return -1;

  int nKeep_i = clash_build_keep_indices(atoms, nAtoms, idx_keep, nAtoms);
  if (nKeep_i <= 0){
    free(idx_keep);
    return -1;
  }
  size_t nKeep = (size_t)nKeep_i;

  for (int it=0; it<=O.max_iter; it++){
    if (place_fn(pts, K, place_user, out_coords) != 0){
      free(idx_keep);
      return -1;
    }

    ClashReport cr;
    if (clash_min_dist_windowed_sortgrid(out_coords, K, nAtoms,
                                         idx_keep, nKeep,
                                         O.W, O.cutoff,
                                         &cr) != 0){
      free(idx_keep);
      return -1;
    }

    R.last_dmin = cr.dmin;

    if (!(cr.dmin < O.cutoff)){
      R.converged = 1;
      R.iters = it;
      break;
    }

    if (it == O.max_iter){
      R.converged = 0;
      R.iters = it;
      break;
    }

    /* grow curve uniformly and retry */
    curve_scale_pts(pts, K, O.grow);
    R.final_scale_mult *= O.grow;
  }

  free(idx_keep);
  if (rep) *rep = R;
  return 0;
}
