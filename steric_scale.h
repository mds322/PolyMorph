#ifndef STERIC_SCALE_H
#define STERIC_SCALE_H

#include "vecops.h"
#include "pdb.h"
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  int W;              /* default 6 */
  double cutoff;      /* default 1.2 */
  double grow;        /* default 1.03 */
  int max_iter;       /* default 20 */
} StericScaleOpt;

typedef struct {
  int iters;
  double final_scale_mult; /* multiplicative factor applied by this module */
  double last_dmin;        /* Angstrom */
  int converged;           /* 1 if no clash at end */
} StericScaleReport;

/* Callback to (re)build the polymer coordinates in out[] given pts[].
   Should return 0 on success. */
typedef int (*StericPlaceFn)(const Vec3 *pts, int K,
                             void *place_user,
                             Vec3 *out_coords);

/* Run steric feedback scaling on pts[] (in-place scaling).
   - atoms/nAtoms used only to derive a heavy-atom index list by name heuristic
   - place_fn called each iteration to compute out_coords
   Returns 0 on success, -1 on error. */
int steric_scale_apply(Vec3 *pts, int K,
                       const Atom *atoms, size_t nAtoms,
                       StericPlaceFn place_fn, void *place_user,
                       Vec3 *out_coords,
                       const StericScaleOpt *opt,
                       StericScaleReport *rep);

#ifdef __cplusplus
}
#endif
#endif
