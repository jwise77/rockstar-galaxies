#ifndef _BOUNDS_H_
#define _BOUNDS_H_

int _check_bounds(float *pos_i, float *pos_f, float *bounds);
int _check_bounds_raw(float *pos_i, float *bounds);
int bounds_overlap(float *b1, float *b2, float *b3, double overlap);
void wrap_into_box(float *pos);
void render_ccache(int64_t *ccache, int64_t ccache_size, float *bounds, float *ccache_bounds, int64_t val);
int64_t check_ccache(int64_t *ccache, int64_t ccache_size, float *pos, float *ccache_bounds);

#endif /* _BOUNDS_H_ */
