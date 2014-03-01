#ifndef _MERGER_H_
#define _MERGER_H_

#include <stdint.h>
#include "halo.h"

struct particle_halo {
  int64_t pid, hid;
};

extern struct halo *halos1, *halos2;
extern struct binary_output_header head1, head2;
extern int64_t *part1, *part2;

void clear_merger_tree(void);
void init_descendants(void);
void connect_particle_ids_to_halo_ids(void);
void calculate_descendants(void);

#endif /* _MERGER_H_ */
