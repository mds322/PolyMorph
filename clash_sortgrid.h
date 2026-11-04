#ifndef CLASH_SORTGRID_H
#define CLASH_SORTGRID_H

#include "vecops.h"
#include "pdb.h"
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  int W;
  double cutoff;
  double dmin;          /* Angstrom, HUGE_VAL if none */
  size_t pairs_tested;
  int early_clash;
} ClashReport;

/* Build indices of atoms to include (heavy atoms by name heuristic). */
int clash_build_keep_indices(const Atom *atoms, size_t nAtoms,
                             int *idx_out, size_t max_out);

/* No-hash grid: sort window atoms by packed cell key and query neighbor cells by binary search.
   out: coords (K*nAtoms) in out[r*nAtoms + i]
   idx_keep/nKeep: indices within 0..nAtoms-1 to include
   W: number of previous repeats to consider
   cutoff: clash cutoff
   Returns 0 on success, -1 on error. */
int clash_min_dist_windowed_sortgrid(const Vec3 *out, int K, size_t nAtoms,
                                     const int *idx_keep, size_t nKeep,
                                     int W, double cutoff,
                                     ClashReport *rep);

#ifdef __cplusplus
}
#endif
#endif
