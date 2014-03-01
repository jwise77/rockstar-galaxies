#ifndef _INTERLEAVING_H_
#define _INTERLEAVING_H_
#include "particle.h"
#include "groupies.h"
#include "rockstar.h"

struct bparticle {
  int64_t id;
  float pos[6];
  int64_t bgid;
  int64_t chunk;
};

struct bgroup {
  int64_t id, chunk, num_p, tagged;
  int64_t next, head;
};

extern int64_t num_new_bp;
extern struct bgroup *bg;
struct bgroup *final_bg;
extern int64_t num_bg_sets;
extern int64_t *bg_set_sizes;
extern int64_t *bg_set_indices;
extern int64_t our_chunk;

void set_bp_chunk(int64_t chunk);
void clear_bp_data(void);
int64_t find_bgroup(struct bparticle *tbp);
void build_bgroup_links(void);
struct bgroup *find_bgroup_from_id(int64_t id, int64_t chunk);
void find_bgroup_sets(int64_t chunk, int64_t *num_sets, int64_t **set_sizes, struct bgroup **groups, int64_t *total_groups);
int64_t calc_next_bgroup_chunk(void);
int64_t prune_setlist(void);
void bgroups_to_setlist(void);
void clear_bg_data(void);
void clear_final_bg_data(void);
void sort_out_halos_for_chunk(int64_t chunk, float *bounds, struct workunit_info *w, struct fof **c_fofs, struct halo **c_halos, struct extra_halo_info **c_ei, struct particle **c_p, struct fof *fofs);
void check_bgroup_sanity(int64_t num_sets, int64_t *set_sizes, struct bgroup *groups);


/* Method
Main analysis process:
*1) Build Smallfofs
*2) Find boundary particles (outside box)
*3) Tag fofs with boundary particles
*4) Build FOFs
*5) Share boundary particles
*6) From known boundary particles in current chunk, check which boundary particles in other chunks match up.
*7) Find unique gids in other chunks to request, associated with each boundary fof
*8) Sort fofs by size.

Client analysis process:
1) Request fof.
2) If boundary fof, request gid lists in other chunks corresponding to the collected gids, lock requested gids. If another client with lower chunk number has asked for a lock, then refuse.  Otherwise, allow.
3) Repeat with the new gids.
4) Ask for lock on particles from each chunk.  
5) Transfer locked particles
6) Run analysis
7) Sort particles into regular particles and boundary particles


8) Boundary particles stored in new location for each main analysis task.

BGC2:
1) Build table of requests for if halo is close to boundary.
2) Exchange requests.
3) Calculate masses.

Temporal Information:
1) Only halos cached; all particle IDs loaded on demand.

MT: No change.
*/


#endif /*_INTERLEAVING_H_*/
