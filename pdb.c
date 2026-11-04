/*PDB Handler*/
#include "pdb.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ---- fixed-column helpers (PDB is 1-based columns) ---- */
static void slice(char *dst, size_t dstsz, const char *src, int col1, int col2){
  int len = col2 - col1 + 1;
  if (len < 0) len = 0;
  if ((size_t)len >= dstsz) len = (int)dstsz - 1;

  int srclen = (int)strlen(src);
  int start = col1 - 1;
  if (start < 0) start = 0;

  if (start >= srclen){
    dst[0] = '\0';
    return;
  }

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

int pdb_read(const char *fname, Atom *atoms, size_t maxAtoms, size_t *outN){
  FILE *fp = fopen(fname, "r");
  if (!fp){ perror("fopen"); return -1; }

  size_t nAtoms = 0;
  char line[PDB_MAX_LINE];

  while (fgets(line, (int)sizeof(line), fp)){
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
    /* normalize record field to either "ATOM " or "HETATM" */
    if (strncmp(a->rec, "HETATM", 6) == 0) strcpy(a->rec, "HETATM");
    else strcpy(a->rec, "ATOM ");

    slice(buf, sizeof(buf), line, 7, 11); trim(buf); a->serial = (*buf)? atoi(buf) : 0;

    slice(a->name, sizeof(a->name), line, 13, 16); /* keep spacing but NUL-term */
    slice(buf, sizeof(buf), line, 17, 17); a->altLoc = buf[0] ? buf[0] : ' ';

    slice(a->resName, sizeof(a->resName), line, 18, 21); trim(a->resName);

    slice(buf, sizeof(buf), line, 22, 22); a->chainID = buf[0] ? buf[0] : ' ';

    slice(buf, sizeof(buf), line, 23, 26); trim(buf); a->resSeq = (*buf)? atoi(buf) : 0;

    slice(buf, sizeof(buf), line, 27, 27); a->iCode = buf[0] ? buf[0] : ' ';

    slice(buf, sizeof(buf), line, 31, 38); trim(buf); a->coord.x = (*buf)? atof(buf) : 0.0;
    slice(buf, sizeof(buf), line, 39, 46); trim(buf); a->coord.y = (*buf)? atof(buf) : 0.0;
    slice(buf, sizeof(buf), line, 47, 54); trim(buf); a->coord.z = (*buf)? atof(buf) : 0.0;

    /* occupancy/tempFactor are optional-ish, but usually present */
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

    nAtoms++;
  }

  fclose(fp);
  *outN = nAtoms;
  return 0;
}

/* ---------------- stats helpers ---------------- */
static int cmp_dbl(const void *a, const void *b){
  double da = *(const double*)a, db = *(const double*)b;
  return (da < db) ? -1 : (da > db) ? 1 : 0;
}

double pdb_median_consecutive_distance(const Vec3 *p, size_t n){
  if (n < 2) return 0.0;

  size_t m = n - 1;
  double *d = (double*)malloc(m * sizeof(double));
  if (!d) return 0.0;

  for (size_t i=0;i<m;i++) d[i] = vec3_dist(p[i+1], p[i]);

  qsort(d, m, sizeof(double), cmp_dbl);

  double med = (m % 2) ? d[m/2] : 0.5*(d[m/2 - 1] + d[m/2]);
  free(d);
  return med;
}

/* ---------------- write PDB ---------------- */
static void pad_atom_name(char out4[5], const char in4[5]){
  size_t n = strlen(in4);
  for (int i=0;i<4;i++) out4[i] = (i<(int)n) ? in4[i] : ' ';
  out4[4] = '\0';
}

/* EXACTLY your existing formatted writer, factored for reuse */
static void pdb_fprint_atom_like_before(FILE *fp, const Atom *a0, Vec3 p){
  char name4[5];
  pad_atom_name(name4, a0->name);

  fprintf(fp, "%-6s%5d %-4s%1c-%4s%1c%4d%1c %8.3f%8.3f%8.3f",
    a0->rec,
    a0->serial,
    name4,
    a0->altLoc ? a0->altLoc : ' ',
    a0->resName,
    a0->chainID ? a0->chainID : ' ',
    a0->resSeq,
    a0->iCode ? a0->iCode : ' ',
    p.x, p.y, p.z);

  if (a0->has_occ_bfac) fprintf(fp, "%6.2f%6.2f", a0->occ, a0->bfac);
  else                 fprintf(fp, "%6s%6s", "1.00", "0.00");

  if (a0->has_elem || a0->has_charge){
    fprintf(fp, "%10s", ""); /* cols 67-76 */
    fprintf(fp, "%2s%2s",
      a0->has_elem ? a0->element : " ",
      a0->has_charge ? a0->charge : " ");
  }
  fputc('\n', fp);
}

int pdb_write_atoms(const char *fname, const Atom *meta_all, const Vec3 *coords_all, size_t n){
  FILE *fp = fopen(fname, "w");
  if (!fp){ perror("fopen"); return -1; }

  for (size_t i=0;i<n;i++){
    pdb_fprint_atom_like_before(fp, &meta_all[i], coords_all[i]);
  }

  fprintf(fp, "END\n");
  fclose(fp);
  return 0;
}

static int same_residue_id(const Atom *a, const Atom *b){
  return (a->chainID == b->chainID) && (a->resSeq == b->resSeq) && (a->iCode == b->iCode);
}

int pdb_count_residues_ordered(const Atom *atoms, size_t nAtoms){
  if (!atoms || nAtoms == 0) return 0;
  int nres = 1;
  for (size_t i=1;i<nAtoms;i++){
    if (!same_residue_id(&atoms[i], &atoms[i-1])) nres++;
  }
  return nres;
}

int pdb_build_residue_index_map(const Atom *atoms, size_t nAtoms, int *out_res_index){
  if (!atoms || nAtoms == 0 || !out_res_index) return -1;

  int cur = 0;
  out_res_index[0] = 0;

  for (size_t i=1;i<nAtoms;i++){
    if (!same_residue_id(&atoms[i], &atoms[i-1])) cur++;
    out_res_index[i] = cur;
  }
  return cur + 1;
}

/* ---- CONECT: link repeats only (SEL2 in unit r -> SEL1 in unit r+1) ---- */
int pdb_write_conect_repeats(FILE *fp,
                             const Atom *meta_all,
                             size_t unitAtoms,
                             int K,
                             int sel1_index0,
                             int sel2_index0)
{
  if (!fp || !meta_all) return -1;
  if (unitAtoms == 0 || K < 2) return 0;
  if (sel1_index0 < 0 || sel2_index0 < 0) return -1;
  if ((size_t)sel1_index0 >= unitAtoms) return -1;
  if ((size_t)sel2_index0 >= unitAtoms) return -1;

  for (int r = 0; r < K-1; r++){
    size_t baseA = (size_t)r * unitAtoms;
    size_t baseB = (size_t)(r+1) * unitAtoms;

    int a = meta_all[baseA + (size_t)sel2_index0].serial;
    int b = meta_all[baseB + (size_t)sel1_index0].serial;

    if (a > 0 && b > 0){
      fprintf(fp, "CONECT%5d%5d\n", a, b);
      fprintf(fp, "CONECT%5d%5d\n", b, a);
    }
  }
  return 0;
}

int pdb_write_atoms_with_conect_repeats(const char *fname,
                                       const Atom *meta_all,
                                       const Vec3 *coords_all,
                                       size_t totalN,
                                       size_t unitAtoms,
                                       int K,
                                       int sel1_index0,
                                       int sel2_index0)
{
  FILE *fp = fopen(fname, "w");
  if (!fp){ perror("fopen"); return -1; }

  for (size_t i=0;i<totalN;i++){
    pdb_fprint_atom_like_before(fp, &meta_all[i], coords_all[i]);
  }

  if (unitAtoms > 0 && totalN == unitAtoms * (size_t)K){
    (void)pdb_write_conect_repeats(fp, meta_all, unitAtoms, K, sel1_index0, sel2_index0);
  }

  fprintf(fp, "END\n");
  fclose(fp);
  return 0;
}

