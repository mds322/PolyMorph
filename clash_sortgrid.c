#include "clash_sortgrid.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef struct {
  uint64_t key;
  Vec3 pos;
} CellAtom;

static void trim_name(char *dst, size_t dstsz, const char name4[5]){
  size_t j=0;
  for (size_t i=0;i<4 && j+1<dstsz;i++){
    char c = name4[i];
    if (c!=' ' && c!='\t' && c!='\r' && c!='\n') dst[j++] = c;
  }
  dst[j]='\0';
}

static int is_hydrogen_name(const Atom *a){
  char t[8];
  trim_name(t, sizeof(t), a->name);
  if (!t[0]) return 0;
  if (t[0] == 'H') return 1;
  if ((t[0] >= '0' && t[0] <= '9') && t[1] == 'H') return 1;
  return 0;
}

int clash_build_keep_indices(const Atom *atoms, size_t nAtoms,
                             int *idx_out, size_t max_out)
{
  if (!atoms || !idx_out) return -1;

  size_t n = 0;
  for (size_t i=0;i<nAtoms;i++){
    if (is_hydrogen_name(&atoms[i])) continue;
    if (n >= max_out) break;
    idx_out[n++] = (int)i;
  }
  if (n == 0){
    for (size_t i=0;i<nAtoms && i<max_out;i++) idx_out[n++] = (int)i;
  }
  return (int)n;
}

static int cmp_cellatom(const void *A, const void *B){
  const CellAtom *a = (const CellAtom*)A;
  const CellAtom *b = (const CellAtom*)B;
  if (a->key < b->key) return -1;
  if (a->key > b->key) return 1;
  return 0;
}

static size_t lower_bound_key(const CellAtom *arr, size_t n, uint64_t key){
  size_t lo=0, hi=n;
  while (lo<hi){
    size_t mid = lo + (hi-lo)/2;
    if (arr[mid].key < key) lo = mid+1;
    else hi = mid;
  }
  return lo;
}

static size_t upper_bound_key(const CellAtom *arr, size_t n, uint64_t key){
  size_t lo=0, hi=n;
  while (lo<hi){
    size_t mid = lo + (hi-lo)/2;
    if (arr[mid].key <= key) lo = mid+1;
    else hi = mid;
  }
  return lo;
}

static int floor_div(double x, double cell){
  return (int)floor(x / cell);
}

/* Pack local (non-negative) cell coords into 64-bit.
   We assume each axis span fits within 21 bits (< 2^21).
   key = ix | (iy<<21) | (iz<<42) */
static uint64_t pack3_21(unsigned ix, unsigned iy, unsigned iz){
  return (uint64_t)ix | ((uint64_t)iy << 21) | ((uint64_t)iz << 42);
}

int clash_min_dist_windowed_sortgrid(const Vec3 *out, int K, size_t nAtoms,
                                     const int *idx_keep, size_t nKeep,
                                     int W, double cutoff,
                                     ClashReport *rep)
{
  if (!out || K < 2 || nAtoms == 0) return -1;
  if (!idx_keep || nKeep == 0) return -1;
  if (!(cutoff > 0.0)) return -1;
  if (W < 1) W = 1;
  if (W > K-1) W = K-1;

  ClashReport R;
  R.W = W;
  R.cutoff = cutoff;
  R.dmin = HUGE_VAL;
  R.pairs_tested = 0;
  R.early_clash = 0;

  const double cell = cutoff;
  const double cut2 = cutoff * cutoff;
  double dmin2 = HUGE_VAL;

  const size_t Nwin_max = (size_t)W * nKeep;
  CellAtom *win = (CellAtom*)malloc(Nwin_max * sizeof(CellAtom));
  if (!win) return -1;

  for (int r=0; r<K; r++){
    int r0 = r - W;
    if (r0 < 0) r0 = 0;

    /* Determine local cell-coordinate bounds for this window */
    int ix_min=0, iy_min=0, iz_min=0;
    int ix_max=0, iy_max=0, iz_max=0;
    int first = 1;

    for (int q=r0; q<r; q++){
      const Vec3 *B = out + (size_t)q * nAtoms;
      for (size_t jj=0; jj<nKeep; jj++){
        int j = idx_keep[jj];
        Vec3 pj = B[(size_t)j];
        int ix = floor_div(pj.x, cell);
        int iy = floor_div(pj.y, cell);
        int iz = floor_div(pj.z, cell);
        if (first){
          ix_min=ix_max=ix; iy_min=iy_max=iy; iz_min=iz_max=iz;
          first = 0;
        } else {
          if (ix < ix_min) ix_min = ix; if (ix > ix_max) ix_max = ix;
          if (iy < iy_min) iy_min = iy; if (iy > iy_max) iy_max = iy;
          if (iz < iz_min) iz_min = iz; if (iz > iz_max) iz_max = iz;
        }
      }
    }
    if (first) continue; /* nothing inserted */

    /* expand bounds by 1 so neighbor queries remain in-range */
    ix_min -= 1; iy_min -= 1; iz_min -= 1;
    ix_max += 1; iy_max += 1; iz_max += 1;

    unsigned span_x = (unsigned)(ix_max - ix_min + 1);
    unsigned span_y = (unsigned)(iy_max - iy_min + 1);
    unsigned span_z = (unsigned)(iz_max - iz_min + 1);

    /* sanity: require spans fit within 21 bits */
    if (span_x >= (1u<<21) || span_y >= (1u<<21) || span_z >= (1u<<21)){
      /* window is too spread out; fall back to increasing cutoff/cell or reduce W.
         For now, treat as error. */
      free(win);
      return -1;
    }

    /* build window list with packed keys */
    size_t nwin = 0;
    for (int q=r0; q<r; q++){
      const Vec3 *B = out + (size_t)q * nAtoms;
      for (size_t jj=0; jj<nKeep; jj++){
        int j = idx_keep[jj];
        Vec3 pj = B[(size_t)j];
        int ix = floor_div(pj.x, cell) - ix_min;
        int iy = floor_div(pj.y, cell) - iy_min;
        int iz = floor_div(pj.z, cell) - iz_min;
        win[nwin].key = pack3_21((unsigned)ix, (unsigned)iy, (unsigned)iz);
        win[nwin].pos = pj;
        nwin++;
      }
    }

    qsort(win, nwin, sizeof(win[0]), cmp_cellatom);

    /* query atoms in repeat r */
    const Vec3 *A = out + (size_t)r * nAtoms;
    for (size_t ii=0; ii<nKeep; ii++){
      int i = idx_keep[ii];
      Vec3 pi = A[(size_t)i];

      int qix0 = floor_div(pi.x, cell) - ix_min;
      int qiy0 = floor_div(pi.y, cell) - iy_min;
      int qiz0 = floor_div(pi.z, cell) - iz_min;

      for (int dz=-1; dz<=1; dz++)
      for (int dy=-1; dy<=1; dy++)
      for (int dx=-1; dx<=1; dx++){
        int qix = qix0 + dx;
        int qiy = qiy0 + dy;
        int qiz = qiz0 + dz;
        if (qix < 0 || qiy < 0 || qiz < 0) continue;
        if ((unsigned)qix >= span_x || (unsigned)qiy >= span_y || (unsigned)qiz >= span_z) continue;

        uint64_t key = pack3_21((unsigned)qix, (unsigned)qiy, (unsigned)qiz);
        size_t lo = lower_bound_key(win, nwin, key);
        size_t hi = upper_bound_key(win, nwin, key);

        for (size_t k=lo; k<hi; k++){
          double d2 = vec3_dist2(pi, win[k].pos);
          R.pairs_tested++;
          if (d2 < dmin2) dmin2 = d2;
          if (d2 < cut2){
            R.early_clash = 1;
            R.dmin = sqrt(dmin2);
            if (rep) *rep = R;
            free(win);
            return 0;
          }
        }
      }
    }
  }

  R.dmin = (dmin2 < HUGE_VAL/2) ? sqrt(dmin2) : HUGE_VAL;
  if (rep) *rep = R;
  free(win);
  return 0;
}
