#ifndef READ_TREE_H
#define READ_TREE_H

#include <stdint.h>

#ifndef EXTRA_HALO_INFO
#define EXTRA_HALO_INFO
#endif

struct halo {
  float scale;
  int64_t id, num_prog, phantom, pid, upid, mmp;
  struct halo *desc, *parent, *uparent, *prog, *next_coprog;
  float mvir, orig_mvir, rvir, rs, vrms, vmax_r, scale_of_last_MM,
    vmax, pos[3], vel[3];
  EXTRA_HALO_INFO
};

struct halo_index_key {
  int64_t id;
  int64_t index;
};

struct halo_list {
  struct halo *halos;
  int64_t num_halos;
  float scale;  
  struct halo_index_key *halo_lookup;
};

struct halo_tree {
  struct halo_list *halo_lists;
  int64_t num_lists;
  int64_t *scale_factor_conv;
  int64_t num_scales;
};

extern struct halo_tree halo_tree;
extern struct halo_list all_halos;

struct halo *lookup_halo_in_list(struct halo_list *hl, int64_t id);
struct halo_list *lookup_scale(float scale);
struct halo_list *find_closest_scale(float scale);
void read_tree(char *filename);
void delete_tree();

#endif /* READ_TREE_H */
