#include <math.h>
#include <assert.h>
#include "config_vars.h"

int _check_bounds(float *pos_i, float *pos_f, float *bounds) {
  int64_t i;
  for (i=0; i<3; i++) {
    if (pos_i[i]>bounds[i+3]) {
      if (bounds[i]<0 && PERIODIC) {
        pos_f[i]=pos_i[i]-BOX_SIZE;
        if (pos_f[i]>bounds[i+3] || pos_f[i] < bounds[i]) return 0;
      } else return 0;
    }
    else if (pos_i[i]<bounds[i]) {
      if (bounds[i+3]>BOX_SIZE && PERIODIC) {
        pos_f[i]=pos_i[i]+BOX_SIZE;
        if (pos_f[i]>bounds[i+3] || pos_f[i] < bounds[i]) return 0;
      } else return 0;
    }
    else pos_f[i] = pos_i[i];
  }
  return 1;
}

void wrap_into_box(float *pos) {
  int64_t i;
  if (!PERIODIC || !BOX_SIZE) return;
  for (i=0; i<3; i++) {
    if (pos[i]>BOX_SIZE) pos[i] -= BOX_SIZE;
    else if (pos[i]<0) pos[i] += BOX_SIZE;
  }
}

int _check_bounds_raw(float *pos_i, float *bounds) {
  int64_t i;
  for (i=0; i<3; i++)
    if (pos_i[i]>=bounds[i+3] || pos_i[i]<bounds[i]) return 0;
  return 1;
}

//Note that b2 is extended by "overlap" in all directions
int bounds_overlap(float *b1, float *b2, float *b3, double overlap) {
  int64_t i;
  int wrap = 0, max_wrap;
  float min, max;
  for (i=0; i<3; i++) {
    b3[i] = min = b2[i]-overlap;
    b3[i+3] = max = b2[i+3]+overlap;
    wrap = ((min < 0) && PERIODIC) ? -1 : 0;
    max_wrap = ((max > BOX_SIZE) && PERIODIC) ? 2 : 1;
    for (; wrap<max_wrap; wrap++) {
      if (((b1[i]+wrap*BOX_SIZE) < max) &&
          ((b1[i+3]+wrap*BOX_SIZE) > min)) break;
    }
    if (wrap == max_wrap) return 0;
  }
  return 1;
}
