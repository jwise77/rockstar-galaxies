#ifndef _SUBHALO_METRIC_H_
#define _SUBHALO_METRIC_H_
#include "halo.h"
#include "particle.h"

struct halo_metric {
  float pos[4]; //Last coordinate is radius
  struct halo *target;
};

void build_subtree(struct halo **subs, int64_t num_subs);
struct halo *find_best_halo(struct particle *part, struct halo *best_halo);
struct halo *find_best_parent(struct halo *h, struct halo *biggest_halo);
float calc_particle_dist(struct halo *h, struct particle *part);
float _calc_halo_dist(struct halo *h1, struct halo *h2);
void free_subtree(void);
struct halo_metric **find_children(struct halo *h, struct halo *parent, 
				   float r, int64_t *num_children);

float calc_halo_dist(struct halo *h1, struct halo *h2);
#endif /* _SUBHALO_METRIC_H_ */
