#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "inthash.h"
#include "check_syscalls.h"
#include "halo.h"
#include "io/meta_io.h"

struct halo *halos1 = NULL, *halos2 = NULL;
struct binary_output_header head1, head2;
int64_t *part1 = NULL, *part2 = NULL;
struct inthash *part2_halos = NULL;
int64_t *part1_halos = NULL;

void clear_merger_tree(void) {
  memset(&head1, 0, sizeof(struct binary_output_header));
  memset(&head2, 0, sizeof(struct binary_output_header));
  halos1 = check_realloc(halos1, 0, "Freeing halos.");
  halos2 = check_realloc(halos2, 0, "Freeing halos.");
  part1 = check_realloc(part1, 0, "Freeing particle IDs.");
  part2 = check_realloc(part2, 0, "Freeing particle IDs.");
  part1_halos = check_realloc(part1_halos, 0, "Freeing ID assignments.");
  if (part2_halos) {
    free_inthash(part2_halos);
    part2_halos = NULL;
  }
}

void init_descendants(void) {
  int64_t i;
  for (i=0; i<head1.num_halos; i++) halos1[i].desc = -1;
}

int compare_int64(const void *a, const void *b) {
  int64_t c = *((int64_t *)a);
  int64_t d = *((int64_t *)b);
  if (c < d) return -1;
  if (c > d) return 1;
  return 0;
}

void connect_particle_ids_to_halo_ids(void) {
  int64_t i,j;
  part2_halos = new_inthash();
  if (!head2.num_particles || !head2.num_halos) return;
  ih_prealloc(part2_halos, head2.num_particles);
  for (i=0; i<head2.num_halos; i++) {
    for (j=halos2[i].p_start; j<halos2[i].p_start+halos2[i].num_p; j++)
      ih_setint64(part2_halos, part2[j], halos2[i].id);
  }
  part2 = check_realloc(part2, 0, "Freeing particle IDs.");  
}

void calculate_descendants(void) {
  int64_t i, j, k, p2;
  int64_t max_p = 0, desc, desc_maxp, last_desc;
  if (!halos2) return;

  for (i=0; i<head1.num_halos; i++) {
    if (halos1[i].num_p > max_p) {
      max_p = halos1[i].num_p;
      part1_halos = check_realloc(part1_halos, sizeof(int64_t)*max_p,
				  "Allocating room for ID assignments.");
    }

    k=0;
    desc = desc_maxp = -1;
    for (j=halos1[i].p_start; j<halos1[i].p_start+halos1[i].num_p; j++) {
      p2 = ih_getint64(part2_halos, part1[j]);
      if (p2 != IH_INVALID) {
	part1_halos[k] = p2;
	k++;
      }
    }

    if (!k) continue;

    qsort(part1_halos, k, sizeof(int64_t), compare_int64);
    last_desc = 0;
    for (j=1; j<k; j++) {
      if (part1_halos[j]!=part1_halos[last_desc]) {
	if (j-last_desc > desc_maxp) {
	  desc_maxp = j - last_desc;
	  desc = part1_halos[last_desc];
	}
	last_desc = j;
      }
    }
    if (j - last_desc > desc_maxp) desc = part1_halos[last_desc];
    halos1[i].desc = desc;
  }
}
