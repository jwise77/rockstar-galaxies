#ifndef _IO_BGC2_H_
#define _IO_BGC2_H_

#include <stdint.h>

extern char **bgc2_snapnames;
extern int64_t num_bgc2_snaps;

struct extended_particle {
  int64_t id, hid;
  float pos[6];
  float mass, energy;
  float softening; /* Per-particle softening, not currently used. */
  float metallicity; /* Not currently used. */
  int32_t type;
};

struct sphere_request {
  float cen[3], r;
};

#define BGC2_R (1.1/1e3)

extern struct extended_particle *ep, *ep2;
extern int64_t num_ep2;

void output_bgc2(int64_t id_offset, int64_t snap, int64_t chunk, float *bounds);
void calc_bgc2_parents(int64_t snap);
int64_t check_bgc2_snap(int64_t snap);
void init_extended_particle_tree(void);
void free_extended_particle_tree(void);
struct extended_particle **do_sphere_request(float cen[3], float r, int64_t *num_sp);

#endif /* _IO_BGC2_H_ */
