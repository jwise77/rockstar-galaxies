/* The Rockstar Halo Finder.
   Copyright (C) 2011-2013  Peter Behroozi

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
   If so, it should be in a file called "LICENSE".
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include "config_vars.h"
#include "config.h"
#include "rockstar.h"
#include "io/meta_io.h"
#include "server.h"
#include "client.h"
#include "version.h"

int main(int argc, char **argv)
{
  int64_t i, s, snap=-1, block=-1, did_config = 0;
  char buffer[1024];
  srand(1);
  if (argc < 2) {
    printf("Rockstar Halo Finder with Galaxy Support, Version %s\n", ROCKSTAR_VERSION);
    printf("(C) 2011-2013 Peter Behroozi.  See the LICENSE file for redistribution details.\n");
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

  if (!PARALLEL_IO) {
    for (i=1; i<argc; i++) {
      if (!strcmp("-c", argv[i])) i++;
      else if (!strcmp("-s", argv[i])) {
	snap = argv[i+1] ? atoi(argv[i+1]) : 0;
	i++;
      }
      else if (!strcmp("-b", argv[i])) {
	block = argv[i+1] ? atoi(argv[i+1]) : 0;
	i++;
      }
      else read_particles(argv[i]);
    }

    if (snap >= 0) {
      if (block < 0) {
	for (block=0; block<NUM_BLOCKS; block++) {
	  get_input_filename(buffer, 1024, snap, block);
	  read_particles(buffer);
	}
      }
      else {
	get_input_filename(buffer, 1024, snap, block);
	read_particles(buffer);
      }
    }

    output_config(NULL);
    rockstar(NULL, 0);
    if (block==NUM_BLOCKS) block = 0;
    if (block < 0) block = 0;
    if (snap < 0) snap = 0;
    output_halos(0,snap,block,NULL);
    free_halos();
  }
  else {
    if (snap > -1) {
      STARTING_SNAP = snap;
      SINGLE_SNAP = 1;
    }
    if ((NUM_WRITERS != 1) && PERIODIC) check_num_writers();
    else if (NUM_WRITERS==1 && PERIODIC && PARALLEL_IO) {
      fprintf(stderr, "[Warning] Setting PERIODIC=0 since NUM_WRITERS=1.\n");
      fprintf(stderr, "[Warning] To enable periodic boundary conditions, increase NUM_WRITERS to at least 8.\n");
      PERIODIC = 0;
    }
    s = server();
    if (!s) client(-1);
  }
  return 0;
}
