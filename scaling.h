#ifndef SCALING_H
#define SCALING_H

#include "pdb.h"
#include "selection.h"
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Compute target step distance (Angstrom) from the input PDB.
   dist(SEL2 - SEL1) in input units
   Returns >0 on success, 0 on failure. */
double scaling_target_step_from_hits(const SelectionHit *h1, const SelectionHit *h2);

#ifdef __cplusplus
}
#endif

#endif
