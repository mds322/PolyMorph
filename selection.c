//Selection functions for monomer addition sites
//
#include "selection.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*Extra helper for splits*/
static char *split_once(char *s, char delim, char **rest){
  if (!s) return NULL;
  char *p = s;
  while (*p && *p != delim) p++;
  if (*p == delim){
    *p = '\0';
    *rest = p + 1;
  } else {
    *rest = NULL;
  }
  return s;
}


/* minimal fixed-column helpers (duplicated to avoid changing pdb.c API) */
static void slice(char *dst, size_t dstsz, const char *src, int col1, int col2){
  int len = col2 - col1 + 1;
  if (len < 0) len = 0;
  if ((size_t)len >= dstsz) len = (int)dstsz - 1;

  int srclen = (int)strlen(src);
  int start = col1 - 1;
  if (start < 0) start = 0;
  if (start >= srclen){ dst[0] = '\0'; return; }

  int avail = srclen - start;
  if (len > avail) len = avail;

  memcpy(dst, src + start, (size_t)len);
  dst[len] = '\0';
}

static void rtrim(char *s){
  size_t n = strlen(s);
  while (n>0 && (s[n-1]==' '||s[n-1]=='\t'||s[n-1]=='\r'||s[n-1]=='\n'))
    s[--n]='\0';
}
static void ltrim(char *s){
  size_t i=0, n=strlen(s);
  while (i<n && (s[i]==' '||s[i]=='\t')) i++;
  if (i>0) memmove(s, s+i, n-i+1);
}
static void trim(char *s){ rtrim(s); ltrim(s); }

int pdb_read_with_lines(const char *fname,
                        Atom *atoms, size_t maxAtoms, size_t *outN,
                        int *line1_out)
{
  FILE *fp = fopen(fname, "r");
  if (!fp){ perror("fopen"); return -1; }

  size_t nAtoms = 0;
  int line1 = 0;
  char line[PDB_MAX_LINE];

  while (fgets(line, (int)sizeof(line), fp)){
    line1++;

    if (strncmp(line, "ATOM", 4) != 0 && strncmp(line, "HETATM", 6) != 0)
      continue;

    if (nAtoms >= maxAtoms){
      fprintf(stderr, "Too many atoms (max %zu)\n", maxAtoms);
      fclose(fp);
      return -1;
    }

    Atom *a = &atoms[nAtoms];
    memset(a, 0, sizeof(*a));

    char buf[64];

    slice(a->rec, sizeof(a->rec), line, 1, 6);
    if (strncmp(a->rec, "HETATM", 6) == 0) strcpy(a->rec, "HETATM");
    else strcpy(a->rec, "ATOM ");

    slice(buf, sizeof(buf), line, 7, 11); trim(buf); a->serial = (*buf)? atoi(buf) : 0;

    slice(a->name, sizeof(a->name), line, 13, 16);
    slice(buf, sizeof(buf), line, 17, 17); a->altLoc = buf[0] ? buf[0] : ' ';

    slice(a->resName, sizeof(a->resName), line, 18, 21); trim(a->resName);

    slice(buf, sizeof(buf), line, 22, 22); a->chainID = buf[0] ? buf[0] : ' ';

    slice(buf, sizeof(buf), line, 23, 26); trim(buf); a->resSeq = (*buf)? atoi(buf) : 0;

    slice(buf, sizeof(buf), line, 27, 27); a->iCode = buf[0] ? buf[0] : ' ';

    slice(buf, sizeof(buf), line, 31, 38); trim(buf); a->coord.x = (*buf)? atof(buf) : 0.0;
    slice(buf, sizeof(buf), line, 39, 46); trim(buf); a->coord.y = (*buf)? atof(buf) : 0.0;
    slice(buf, sizeof(buf), line, 47, 54); trim(buf); a->coord.z = (*buf)? atof(buf) : 0.0;

    slice(buf, sizeof(buf), line, 55, 60); trim(buf);
    if (*buf){ a->occ = atof(buf); a->has_occ_bfac = 1; }
    slice(buf, sizeof(buf), line, 61, 66); trim(buf);
    if (*buf){ a->bfac = atof(buf); a->has_occ_bfac = 1; }

    slice(a->element, sizeof(a->element), line, 77, 78);
    trim(a->element);
    a->has_elem = (a->element[0] != '\0');

    slice(a->charge, sizeof(a->charge), line, 79, 80);
    trim(a->charge);
    a->has_charge = (a->charge[0] != '\0');

    if (line1_out) line1_out[nAtoms] = line1;
    nAtoms++;
  }

  fclose(fp);
  *outN = nAtoms;
  return 0;
}

static int atom_name_trim_eq(const Atom *a, const char *want){
  char tmp[8];
  snprintf(tmp, sizeof(tmp), "%s", a->name);
  trim(tmp);
  return strcmp(tmp, want) == 0;
}

static int match_chain_res_atom(const Atom *a, char chainID, int resSeq, char iCode, const char *atomName){
  if (a->chainID != chainID) return 0;
  if (a->resSeq  != resSeq)  return 0;
  if (a->iCode   != iCode)   return 0;
  return atom_name_trim_eq(a, atomName);
}

int selection_resolve_one(const Atom *atoms, size_t nAtoms,
                          const int *line1,
                          const SelectionSpec *spec,
                          SelectionHit *hit)
{
  if (!atoms || !spec || !hit) return -1;
  memset(hit, 0, sizeof(*hit));

  int found_idx = -1;

  switch (spec->kind){
    case SEL_BY_NAME_FIRSTN: {
      int nth = spec->nth;
      if (nth < 1) nth = 1;
      int count = 0;
      for (size_t i=0;i<nAtoms;i++){
        if (atom_name_trim_eq(&atoms[i], spec->atomName)){
          count++;
          if (count == nth){ found_idx = (int)i; break; }
        }
      }
    } break;

    case SEL_BY_SERIAL: {
      for (size_t i=0;i<nAtoms;i++){
        if (atoms[i].serial == spec->serial){ found_idx = (int)i; break; }
      }
    } break;

    case SEL_BY_INDEX_0: {
      if (spec->index0 >= 0 && (size_t)spec->index0 < nAtoms) found_idx = spec->index0;
    } break;

    case SEL_BY_LINE_1: {
      if (!line1) return -1;
      for (size_t i=0;i<nAtoms;i++){
        if (line1[i] == spec->line1){ found_idx = (int)i; break; }
      }
    } break;

    case SEL_BY_CHAIN_RESSEQ_ICODE_ATOM: {
      char chainID = spec->chainID ? spec->chainID : ' ';
      char iCode   = spec->iCode   ? spec->iCode   : ' ';
      for (size_t i=0;i<nAtoms;i++){
        if (match_chain_res_atom(&atoms[i], chainID, spec->resSeq, iCode, spec->atomName2)){
          found_idx = (int)i;
          break;
        }
      }
    } break;

    default:
      return -1;
  }

  if (found_idx < 0) return -1;

  hit->ok = 1;
  hit->index0 = found_idx;
  hit->line1  = (line1 ? line1[found_idx] : 0);
  hit->atom   = atoms[found_idx];
  return 0;
}

int selection_resolve_pair(const Atom *atoms, size_t nAtoms,
                           const int *line1,
                           const SelectionSpec *s1,
                           const SelectionSpec *s2,
                           SelectionHit *h1,
                           SelectionHit *h2)
{
  if (selection_resolve_one(atoms, nAtoms, line1, s1, h1) != 0) return -1;
  if (selection_resolve_one(atoms, nAtoms, line1, s2, h2) != 0) return -1;
  return 0;
}

SelectionSpec selection_default_O4_1(void){
  SelectionSpec s;
  memset(&s, 0, sizeof(s));
  s.kind = SEL_BY_NAME_FIRSTN;
  strcpy(s.atomName, "O4");
  s.nth = 1;
  return s;
}

SelectionSpec selection_default_O4_2(void){
  SelectionSpec s;
  memset(&s, 0, sizeof(s));
  s.kind = SEL_BY_NAME_FIRSTN;
  strcpy(s.atomName, "O4");
  s.nth = 2;
  return s;
}

/* -------- CLI parsing helpers -------- */

static int parse_int(const char *s, int *out){
  if (!s || !*s) return -1;
  char *end = NULL;
  long v = strtol(s, &end, 10);
  if (!end || *end != '\0') return -1;
  if (v < -2147483647L || v > 2147483647L) return -1;
  *out = (int)v;
  return 0;
}

int selection_parse(const char *arg, SelectionSpec *out){
  if (!arg || !out) return -1;
  memset(out, 0, sizeof(*out));

  if (strncmp(arg, "name:", 5) == 0){
    const char *p = arg + 5;
    if (!*p) return -1;

    /* copy atom name up to '#' or end */
    char name[8];
    size_t j = 0;
    while (*p && *p != '#' && j+1 < sizeof(name)) name[j++] = *p++;
    name[j] = '\0';
    if (name[0] == '\0') return -1;

    int nth = 1;
    if (*p == '#'){
      p++;
      if (parse_int(p, &nth) != 0 || nth < 1) return -1;
    }

    out->kind = SEL_BY_NAME_FIRSTN;
    strncpy(out->atomName, name, sizeof(out->atomName)-1);
    out->atomName[sizeof(out->atomName)-1] = '\0';
    out->nth = nth;
    return 0;
  }

  if (strncmp(arg, "serial:", 7) == 0){
    int v;
    if (parse_int(arg+7, &v) != 0) return -1;
    out->kind = SEL_BY_SERIAL;
    out->serial = v;
    return 0;
  }

  if (strncmp(arg, "index:", 6) == 0){
    int v;
    if (parse_int(arg+6, &v) != 0 || v < 0) return -1;
    out->kind = SEL_BY_INDEX_0;
    out->index0 = v;
    return 0;
  }

  if (strncmp(arg, "line:", 5) == 0){
    int v;
    if (parse_int(arg+5, &v) != 0 || v < 1) return -1;
    out->kind = SEL_BY_LINE_1;
    out->line1 = v;
    return 0;
  }

  if (strncmp(arg, "res:", 4) == 0){
    /* res:CHAIN:RESSEQ:ICODE:ATOM
       CHAIN: 1 char; '_' => space
       ICODE: 1 char; '_' => space
     */
    char buf[128];
    strncpy(buf, arg+4, sizeof(buf)-1);
    buf[sizeof(buf)-1] = '\0';

    char *rest = buf;
    char *tok_chain = split_once(rest, ':', &rest);
    char *tok_resseq = rest ? split_once(NULL, ':', &rest) : NULL;
    char *tok_icode  = rest ? split_once(NULL, ':', &rest) : NULL;
    char *tok_atom   = rest ? rest : NULL;
    if (!tok_chain || !tok_resseq || !tok_icode || !tok_atom) return -1;

    int resSeq;
    if (parse_int(tok_resseq, &resSeq) != 0) return -1;

    char chainID = tok_chain[0] ? tok_chain[0] : ' ';
    char iCode   = tok_icode[0] ? tok_icode[0] : ' ';
    if (chainID == '_') chainID = ' ';
    if (iCode   == '_') iCode   = ' ';

    out->kind = SEL_BY_CHAIN_RESSEQ_ICODE_ATOM;
    out->chainID = chainID;
    out->resSeq  = resSeq;
    out->iCode   = iCode;
    strncpy(out->atomName2, tok_atom, sizeof(out->atomName2)-1);
    out->atomName2[sizeof(out->atomName2)-1] = '\0';
    return 0;
  }

  return -1;
}

void selection_print_usage(FILE *fp){
  fprintf(fp,
    "Selection formats:\n"
    "  name:ATOM[#N]            e.g. name:O4#1, name:O4#2\n"
    "  serial:INT               e.g. serial:123\n"
    "  index:INT                0-based atom index in parsed ATOM/HETATM list\n"
    "  line:INT                 1-based input file line number (ATOM/HETATM only)\n"
    "  res:CH:RESSEQ:IC:ATOM    e.g. res:A:12:_:O4  (_ means space)\n");
}

void selection_spec_print(FILE *fp, const SelectionSpec *s){
  if (!s){ fprintf(fp, "(null)\n"); return; }
  switch (s->kind){
    case SEL_BY_NAME_FIRSTN:
      fprintf(fp, "SEL_BY_NAME_FIRSTN name=%s nth=%d\n", s->atomName, s->nth);
      break;
    case SEL_BY_SERIAL:
      fprintf(fp, "SEL_BY_SERIAL serial=%d\n", s->serial);
      break;
    case SEL_BY_INDEX_0:
      fprintf(fp, "SEL_BY_INDEX_0 index0=%d\n", s->index0);
      break;
    case SEL_BY_LINE_1:
      fprintf(fp, "SEL_BY_LINE_1 line1=%d\n", s->line1);
      break;
    case SEL_BY_CHAIN_RESSEQ_ICODE_ATOM:
      fprintf(fp, "SEL_BY_CHAIN_RESSEQ_ICODE_ATOM chain=%c resSeq=%d iCode=%c atom=%s\n",
              s->chainID ? s->chainID : ' ',
              s->resSeq,
              s->iCode ? s->iCode : ' ',
              s->atomName2);
      break;
    default:
      fprintf(fp, "SEL kind=%d (unknown)\n", (int)s->kind);
  }
}
