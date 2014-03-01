#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <string.h>
#include <assert.h>
#include "../config_vars.h"
#include "../config.h"
#include "../rockstar.h"
#include "../io/meta_io.h"
#include "../io/io_bgc2.h"
#include "../check_syscalls.h"
#include "../io/bgc2.h"

extern GROUP_DATA_RMPVMAX *gd;
void load_bgc2_groups(char *filename, struct bgc2_header *hdr,
		      GROUP_DATA_RMPVMAX **groups, int64_t *num_groups);

int main(int argc, char **argv)
{
  int64_t i, snap=-1, did_config = 0;
  char buffer[1024];
  int64_t num_groups = 0;
  struct bgc2_header *hdrs = NULL;

  if (argc < 2) {
    printf("Usage: %s [-c config] [-s snapnum] [-b blocknum]\n", argv[0]);
    exit(1);
  }

  for (i=1; i<argc-1; i++) {
    if (!strcmp("-c", argv[i])) { do_config(argv[i+1]); i++; did_config=1; }
    if (!strcmp("-s", argv[i])) { snap = atoi(argv[i+1]); i++; }
  }
  if (!did_config) do_config(NULL);
  if (strlen(SNAPSHOT_NAMES)) 
    read_input_names(SNAPSHOT_NAMES, &snapnames, &NUM_SNAPS);
  if (strlen(BLOCK_NAMES))
    read_input_names(BLOCK_NAMES, &blocknames, &NUM_BLOCKS);



  hdrs = check_realloc(hdrs, BGC2_HEADER_SIZE*NUM_WRITERS,
		       "Allocating BGC2 headers.");
  for (i=0; i<NUM_WRITERS; i++) {
    get_output_filename(buffer, 1024, snap, i, "bgc2");
    load_bgc2_groups(buffer, hdrs + i, &gd, &num_groups);
  }
  printf("#ID DescID M%s Vmax Vrms R%s Rs Np X Y Z VX VY VZ Parent_ID\n",
	 MASS_DEFINITION, MASS_DEFINITION);
  for (i=0; i<num_groups; i++) {
    printf("%"PRId64" -1 %.5g %.3f 0 %.3f 0 %"PRId64" %f %f %f %f %f %f %"PRId64"\n", gd[i].id, gd[i].mass, gd[i].vmax, gd[i].radius*1.0e3, gd[i].npart, gd[i].pos[0],  gd[i].pos[1],  gd[i].pos[2], gd[i].vel[0], gd[i].vel[1], gd[i].vel[2], gd[i].parent_id);
  }
  return 0;
}
