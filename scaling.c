#include "scaling.h"
#include <math.h>

double scaling_target_step_from_hits(const SelectionHit *h1, const SelectionHit *h2){
  if (!h1 || !h2) return 0.0;
  if (!h1->ok || !h2->ok) return 0.0;

  double dx = h2->atom.coord.x - h1->atom.coord.x;
  double dy = h2->atom.coord.y - h1->atom.coord.y;
  double dz = h2->atom.coord.z - h1->atom.coord.z;
  double d = sqrt(dx*dx + dy*dy + dz*dz);
  return (d > 0.0) ? d : 0.0;
}
