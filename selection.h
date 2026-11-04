/*
 * selection.h
 * Selects pairs to define monomer addition sites
*/

#ifndef SELECTION_H
#define SELECTION_H

#include "pdb.h"
#include <stddef.h>
#include <stdio.h>   /* FILE */

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
  SEL_BY_NAME_FIRSTN = 1,            /* first N atoms matching trimmed name (e.g. O4#2) */
  SEL_BY_SERIAL      = 2,            /* PDB atom serial field */
  SEL_BY_INDEX_0     = 3,            /* 0-based index in parsed ATOM/HETATM array */
  SEL_BY_LINE_1      = 4,            /* 1-based input file line number (ATOM/HETATM only) */
  SEL_BY_CHAIN_RESSEQ_ICODE_ATOM = 5 /* chainID + resSeq + iCode + atom name (trimmed) */
} SelectionKind;

typedef struct {
  SelectionKind kind;

  /* SEL_BY_NAME_FIRSTN */
  char atomName[8];      /* trimmed match, e.g. "O4" */
  int  nth;              /* 1 => first match, 2 => second match, ... */

  /* SEL_BY_SERIAL */
  int serial;

  /* SEL_BY_INDEX_0 */
  int index0;

  /* SEL_BY_LINE_1 */
  int line1;

  /* SEL_BY_CHAIN_RESSEQ_ICODE_ATOM */
  char chainID;          /* ' ' allowed */
  int  resSeq;
  char iCode;            /* ' ' allowed */
  char atomName2[8];     /* trimmed */
} SelectionSpec;

typedef struct {
  int   ok;              /* 1 if found */
  int   index0;          /* index in atoms[] */
  int   line1;           /* original ATOM/HETATM line number (1-based); 0 if unknown */
  Atom  atom;            /* copy of matched atom record */
} SelectionHit;

/* Read PDB and also capture per-atom ATOM/HETATM line numbers.
   line1_out must be maxAtoms long. line1_out[i] set for i<nAtoms. */
int pdb_read_with_lines(const char *fname,
                        Atom *atoms, size_t maxAtoms, size_t *outN,
                        int *line1_out);

/* Resolve a single selection spec to a hit. Returns 0 on success, -1 on failure. */
int selection_resolve_one(const Atom *atoms, size_t nAtoms,
                          const int *line1, /* may be NULL if not available */
                          const SelectionSpec *spec,
                          SelectionHit *hit);

/* Convenience: resolve two selections (e.g., axis definition). */
int selection_resolve_pair(const Atom *atoms, size_t nAtoms,
                           const int *line1,
                           const SelectionSpec *s1,
                           const SelectionSpec *s2,
                           SelectionHit *h1,
                           SelectionHit *h2);

/* Defaults that preserve current behavior: first two atoms named "O4". */
SelectionSpec selection_default_O4_1(void);
SelectionSpec selection_default_O4_2(void);

/* Parse a selection string into SelectionSpec.
   Supported:
     name:ATOM[#N]            e.g. name:O4#2
     serial:INT               e.g. serial:123
     index:INT                0-based
     line:INT                 1-based ATOM/HETATM line number
     res:CH:RESSEQ:IC:ATOM    e.g. res:A:12:_:O4  (_ means space)
   Returns 0 on success, -1 on error. */
int selection_parse(const char *s, SelectionSpec *out);

/* Print selection syntax help. */
void selection_print_usage(FILE *fp);

void selection_spec_print(FILE *fp, const SelectionSpec *s);

#ifdef __cplusplus
}
#endif

#endif

