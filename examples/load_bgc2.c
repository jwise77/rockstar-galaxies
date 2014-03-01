/* Example code for loading BGC2 output files. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <math.h>
#include <assert.h>
#include <sys/stat.h>
#include "../check_syscalls.h"
#include "../io/bgc2.h"
#include "../io/io_util.h"

GROUP_DATA_RMPVMAX *grps = NULL;
PARTICLE_DATA_PVM *parts = NULL;
int64_t num_groups=0, num_parts=0;

void load_bgc2(char *filename, struct bgc2_header *hdr,
	       GROUP_DATA_RMPVMAX **groups, int64_t *num_groups,
	       PARTICLE_DATA_PVM **pdata, int64_t *num_parts);


void load_bgc2(char *filename, struct bgc2_header *hdr,
	       GROUP_DATA_RMPVMAX **groups, int64_t *num_groups,
	       PARTICLE_DATA_PVM **pdata, int64_t *num_parts)
{
  FILE *input;
  int64_t new_group_size, new_part_size;
  int64_t i, p_start;

  assert(sizeof(struct bgc2_header) == BGC2_HEADER_SIZE);
  input = check_fopen(filename, "rb");

  fread_fortran(hdr, BGC2_HEADER_SIZE, 1, input, 0);
  assert(hdr->magic == BGC_MAGIC);
  assert(hdr->version == 2);
  assert(hdr->format_group_data == GDATA_FORMAT_RMPVMAX);

  new_group_size = sizeof(GROUP_DATA_RMPVMAX)*((*num_groups)+hdr->ngroups);
  *groups = check_realloc(*groups, new_group_size, "Allocating groups.");
  fread_fortran((*groups) + (*num_groups), sizeof(GROUP_DATA_RMPVMAX), 
		hdr->ngroups, input, 0);
  
  new_part_size = sizeof(PARTICLE_DATA_PVM)*((*num_parts)+hdr->npart);
  *pdata = check_realloc(*pdata, new_part_size, "Allocating particles");
  p_start = (*num_parts);
  for (i=0; i<hdr->ngroups; i++) {
    fread_fortran((*pdata) + p_start, sizeof(PARTICLE_DATA_PVM), 
		  groups[0][(*num_groups)+i].npart, input, 0);
    p_start += groups[0][(*num_groups)+i].npart;
  }
  *num_groups += hdr->ngroups;
  *num_parts += hdr->npart;
  fclose(input);
}


