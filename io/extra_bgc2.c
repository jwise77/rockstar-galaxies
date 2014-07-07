#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "bgc2.h"
#include "meta_io.h"
#include "io_util.h"
#include "../config_vars.h"
#include "../check_syscalls.h"
#include "../universal_constants.h"
#include "../rockstar.h"
#include "../groupies.h"

extern char **bgc2_snapnames;
extern int64_t num_bgc2_snaps;
extern GROUP_DATA_RMPVMAX *gd;

#define FAST3TREE_TYPE GROUP_DATA_RMPVMAX
#define FAST3TREE_PREFIX EXTRA_BGC2
#include "../fast3tree.c"

#define GROUP_LIST gd
#define RADIUS radius
#define RADIUS_CONVERSION 1.0
#define parent parent_id
#include "../parents.c"
#undef parent

void load_bgc2_groups(char *filename, struct bgc2_header *hdr,
		      GROUP_DATA_RMPVMAX **groups, int64_t *num_groups);

int sort_by_id(const void *a, const void *b) {
  const GROUP_DATA_RMPVMAX *c = a;
  const GROUP_DATA_RMPVMAX *d = b;
  if (c->id < d->id) return -1;
  if (c->id > d->id) return 1;
  return 0;
}

void calc_bgc2_parents(int64_t snap)
{
  char buffer[1024];
  int64_t i, num_groups = 0;
  struct bgc2_header *hdrs = NULL;
  FILE *output;
  int64_t *first_ids = NULL;
  GROUP_DATA_RMPVMAX key, *first_group;

  hdrs = check_realloc(hdrs, BGC2_HEADER_SIZE*NUM_WRITERS,
		       "Allocating BGC2 headers.");
  first_ids = check_realloc(first_ids, sizeof(int64_t)*NUM_WRITERS,
		       "Allocating IDs.");

  for (i=0; i<NUM_WRITERS; i++) {
    get_output_filename(buffer, 1024, snap, i, "bgc2");
    load_bgc2_groups(buffer, hdrs + i, &gd, &num_groups);
    BOX_SIZE = hdrs[i].box_size;
    first_ids[i]=-1;
    if (hdrs[i].ngroups) {
      first_ids[i] = gd[num_groups - hdrs[i].ngroups].id;
    }
    if (!i) hdrs[0].npart_total = hdrs[0].ngroups_total = 0;
    hdrs[0].npart_total += hdrs[i].npart;
    hdrs[0].ngroups_total += hdrs[i].ngroups;
    if (hdrs[i].max_npart > hdrs[0].max_npart_total)
      hdrs[0].max_npart_total = hdrs[i].max_npart;
  }

  for (i=1; i<NUM_WRITERS; i++) {
    hdrs[i].npart_total = hdrs[0].npart_total;
    hdrs[i].ngroups_total = hdrs[0].ngroups_total;
    hdrs[i].max_npart_total = hdrs[0].max_npart_total;
  }

  find_parents(num_groups);
  qsort(gd, num_groups, sizeof(GROUP_DATA_RMPVMAX), sort_by_id);

  for (i=0; i<NUM_WRITERS; i++) {
    if (hdrs[i].ngroups) {
      assert(first_ids[i] != -1);
      key.id = first_ids[i];
      first_group = bsearch(&key, gd, num_groups, sizeof(GROUP_DATA_RMPVMAX), sort_by_id);
      assert(first_group);
    } else {
      first_group = gd;
    }
    get_output_filename(buffer, 1024, snap, i, "bgc2");
    output = fopen(buffer, "r+b");
    fwrite_fortran(hdrs + i, BGC2_HEADER_SIZE, 1, output);
    fwrite_fortran(first_group, sizeof(GROUP_DATA_RMPVMAX), 
		   hdrs[i].ngroups, output);
    fclose(output);
  }

  gd = check_realloc(gd, 0, "Freeing group data.");
  free(first_ids);
  free(hdrs);
}
