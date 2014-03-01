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
#include "load_bgc2.h"


int main(int argc, char **argv) {
  int64_t i,j,k,next_halo_k,halo_sub_k,is_sub;
  struct bgc2_header hdr;

  if (argc < 2) {
    printf("Usage: %s file1.bgc2 ...\n", argv[0]);
    exit(1);
  }

  printf("#X Y Z VX VY VZ Mass ParticleID IsSubstructure HaloID HaloParent\n");
  printf("#Units: positions in Mpc/h (comoving)\n");
  printf("#Units: velocities in km/s (non-comoving)\n");
  printf("#IsSubstructure: 1 if the particle belongs to substructure within the halo; 0 otherwise\n");
  printf("#HaloParent: the ID of a larger halo containing this halo, if any; -1 otherwise\n");
  for (i=1; i<argc; i++) {
    num_groups = num_parts = 0;
    load_bgc2(argv[i], &hdr, &grps, &num_groups, &parts, &num_parts);
    if (i==1) {
      printf("#Box size: %f Mpc/h\n", hdr.box_size);
      printf("#Particle Mass: %f Msun/h\n", hdr.part_mass);
      printf("#h: %f\n", hdr.Hubble0);
      printf("#Redshift: %f\n", hdr.redshift);
    }
    //Do stuff with particles here...
    //See io/bgc2.h for a description of the group and particle structure format
    is_sub=k=j=0;
    next_halo_k = k+grps[j].npart;
    halo_sub_k = k+grps[j].npart_self;
    for (; k<num_parts; k++) {
      if (k==halo_sub_k) is_sub=1;
      if (k==next_halo_k) {
	is_sub=0;
	j++;
	next_halo_k = k+grps[j].npart;
	halo_sub_k = k+grps[j].npart_self;
      }
      printf("%f %f %f %f %f %f %e %"PRId64" %"PRId64" %"PRId64" %"PRId64"\n",
	     parts[k].pos[0], parts[k].pos[1], parts[k].pos[2], 
	     parts[k].vel[0], parts[k].vel[1], parts[k].vel[2], 
	     parts[k].mass, parts[k].part_id, is_sub, grps[j].id, grps[j].parent_id);
    }
  }
  return 0;
}
