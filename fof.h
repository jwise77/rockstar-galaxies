#ifndef FOF_H
#define FOF_H

#include <stdint.h>
#include "particle.h"

struct fof {
  int64_t num_p;
  struct particle *particles;
};

struct smallfof {
  int64_t root;
};

void init_particle_smallfofs(int64_t num_p, struct particle *particles);
void link_particle_to_fof(struct particle *p, int64_t n, struct particle **links);
void link_fof_to_fof(struct particle *p, int64_t n, struct particle **links);
int64_t tag_boundary_particle(struct particle *p);

void build_fullfofs(void);
struct fof *return_fullfofs(int64_t *num_f, int64_t *num_bf);
void copy_fullfofs(struct fof **base, int64_t *num_f, int64_t *num_alloced_f);

void partition_sort_particles(int64_t min, int64_t max,
			      struct particle *particles, int64_t *assignments);

#endif /* FOF_H */
