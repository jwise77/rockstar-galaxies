#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include "../config_vars.h"
#include "../config.h"
#include "../rockstar.h"
#include "../io/meta_io.h"
#include "../io/io_bgc2.h"

int main(int argc, char **argv)
{
  int64_t i, snap=-1, did_config = 0;
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

  calc_bgc2_parents(snap);
  return 0;
}
