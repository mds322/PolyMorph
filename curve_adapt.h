#ifndef CURVE_ADAPT_H
#define CURVE_ADAPT_H

#include "curve.h"
#ifndef M_PI
#define M_PI (3.14159265358979323846264338327950)
#endif /* M_PI */

#ifndef M_E
#define M_E  (2.71828182845904523536028747135266)
#endif /* M_E */


#ifdef __cplusplus
extern "C" {
#endif

/* Max turning angle (radians) between consecutive segment tangents of ws->p[].
   Returns >=0 on success, negative on error. */
double curve_ws_max_turn_angle(const CurveWS *ws);

/* Build curve_ws with adaptive M (caller-owned buffers).
   If max turn angle > theta_limit_rad, request more M by returning +1 and setting ws->M *= 2 (capped at maxM).
   Caller must realloc ws->p/ws->s to new ws->M and call again.
   Returns 0 when quality is acceptable, -1 on error, +1 meaning "need bigger buffers". */
int curve_ws_build_adaptive_turn(CurvePosFn f, void *user,
                                 double t0, double t1,
                                 CurveWS *ws,
                                 int maxM,
                                 double theta_limit_rad);

#ifdef __cplusplus
}
#endif
#endif
