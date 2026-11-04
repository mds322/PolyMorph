/*PDB Handler*/
#ifndef PDB_H
#define PDB_H
#include <stdio.h>
#include <stddef.h>
#include "vecops.h"

#ifdef __cplusplus
extern "C" {
#endif

#define PDB_MAX_ATOMS 999999
#define PDB_MAX_LINE  256

typedef struct {
  char rec[7];      /* "ATOM " or "HETATM" (normalized) */
  int  serial;
  char name[5];     /* cols 13-16 (kept as-is, incl spaces), NUL-term */
  char altLoc;      /* col 17 */
  char resName[5];  /* cols 18-20 trimmed */
  char chainID;     /* col 22 */
  int  resSeq;      /* cols 23-26 */
  char iCode;       /* col 27 */
  Vec3 coord;       /* cols 31-54 */
  double occ;       /* cols 55-60 (optional) */
  double bfac;      /* cols 61-66 (optional) */
  char element[3];  /* cols 77-78 trimmed (optional) */
  char charge[3];   /* cols 79-80 trimmed (optional) */

  int has_occ_bfac;
  int has_elem;
  int has_charge;
} Atom;

/* Read ATOM/HETATM records from a PDB file into atoms[] (up to maxAtoms). */
int  pdb_read(const char *fname, Atom *atoms, size_t maxAtoms, size_t *outN);

/* median of |p[i+1]-p[i]| for i=0..n-2 (returns 0 if n<2) */
double pdb_median_consecutive_distance(const Vec3 *p, size_t n);

/* Write ATOM/HETATM records for n atoms from per-atom metadata + coords. */
int  pdb_write_atoms(const char *fname, const Atom *meta_all, const Vec3 *coords_all, size_t n);

/* Return number of residues in atoms[] when grouped by (chainID,resSeq,iCode)
   in file order. Returns 0 if nAtoms==0. */
int pdb_count_residues_ordered(const Atom *atoms, size_t nAtoms);

/* Build a residue map for atoms[]: out_res_index[i] is 0..(nres-1).
   Residues are defined by (chainID,resSeq,iCode) and assigned in file order.
   Returns nres on success, -1 on error. */
int pdb_build_residue_index_map(const Atom *atoms, size_t nAtoms, int *out_res_index);

/*Add CONECT Records between monomers for visualization*/
int pdb_write_conect_repeats(FILE *fp, const Atom *meta_all, size_t unitAtoms, int K, int sel1_index0, int sel2_index0);

int pdb_write_atoms_with_conect_repeats(const char *fname, const Atom *meta_all, const Vec3 *coords_all, size_t totalN, size_t unitAtoms, int K, int sel1_index0, int sel2_index0);

#ifdef __cplusplus
}
#endif
#endif
