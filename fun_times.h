#ifndef _FUN_TIMES_H_
#define _FUN_TIMES_H_

#include <inttypes.h>

struct previous_halo {
  float pos[6];
  int64_t *particles;
  int64_t num_p, id;
  int64_t chunk;
  float m, r;
};

struct prev_bounds {
  float bounds[6];
};

struct prev_bounds *p_bounds;
extern int64_t prev_snap;

#define MAX_CORE_PARTICLES 10000
#define PREV_HALO_BUFFER_SIZE 100000

void load_previous_halos(int64_t snap, int64_t chunk, float *bounds);
void convert_and_sort_core_particles(struct halo *h, struct particle *hp, float max_r, int64_t *n_core);
float find_previous_mass(struct halo *h, struct particle *hp, int64_t *best_num_p, float max_r);
void clear_prev_files(void);
void reassign_particles_to_parent(struct halo *h, struct particle *hp, int64_t *particle_halos, struct halo *parent_h);

#endif /* _FUN_TIMES_H_ */
