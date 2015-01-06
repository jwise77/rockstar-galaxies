#ifndef GROUPIES_H
#define GROUPIES_H

#include "fof.h"
#include "halo.h"
#include <stdint.h>

#define GROWING_FLAG 1
#define DELETE_FLAG 2
#define POSSIBLE_SWAP_FLAG 4
#define TAGGED_FLAG 8
#define ALWAYS_PRINT_FLAG 16
#define GALAXY_INELIGIBLE_FLAG 32

extern struct halo *halos;
extern int64_t num_halos;
extern struct extra_halo_info *extra_info;

void find_subs(struct fof *f);
void calc_mass_definition(void);
void free_particle_copies(void);
void free_halos(void);
float max_halo_radius(struct halo *h);

//Internal functions
void add_new_halo(void);
void alloc_particle_copies(int64_t total_copies);
void norm_sd(struct fof *f, float thresh, double *axis_x, double *axis_v);
float find_median_r(float *rad, int64_t num_p, float frac);

void find_sd(struct particle *particles, int64_t n, double *sig_x, double *sig_v);
void norm_sd_bary(struct fof *f);

#endif /* GROUPIES_H */
