#ifndef LOAD_BGC2
#define LOAD_BGC2

#include <inttypes.h>
#include "../io/bgc2.h"

extern GROUP_DATA_RMPVMAX *grps;
extern PARTICLE_DATA_PVM *parts;
extern int64_t num_groups, num_parts;

void load_bgc2(char *filename, struct bgc2_header *hdr,
	       GROUP_DATA_RMPVMAX **groups, int64_t *num_groups,
	       PARTICLE_DATA_PVM **pdata, int64_t *num_parts);

#endif /* LOAD_BGC2 */
