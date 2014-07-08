#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <rpc/types.h>
#include <rpc/xdr.h>
#include <string.h>
#include "../check_syscalls.h"
#include "../universal_constants.h"
#include "../config_vars.h"
#include "../config.h"
#include "../particle.h"
#include "io_util.h"
#include "io_nchilada.h"

enum NCDataTypeCode {
  int8 = 1,
  uint8,
  int16,
  uint16,
  int32,
  uint32,
  int64,
  uint64,
  float32,
  float64
};


#define OFFSET(x) ((char *)&(p[0][0].x)-(char*)(p[0]))
#define NCHILADA_RBUFFER 1000000
int64_t NCHILADA_SWAP = 0;
int64_t NCHILADA_ID_BYTES = 4;

int64_t read_nchilada_header(FILE *input, char *filename, struct nchilada_header *nh, int64_t dims) {
  check_fread(&nh->magic, sizeof(int32_t), 1, input);
  if (nh->magic != NCHILADA_MAGIC) {
    NCHILADA_SWAP = 1;
    swap_endian_4byte((int8_t *)&nh->magic);
    if (nh->magic != NCHILADA_MAGIC) {
      fprintf(stderr, "[Error] Failed to read valid NChilada header from file %s!\n", filename);
    }
  }
  fread_mswap(&nh->time, sizeof(double), 1, input, NCHILADA_SWAP);
  fread_mswap(&nh->iHighWord, sizeof(int32_t), 1, input, NCHILADA_SWAP);
  fread_mswap(&nh->nbodies, sizeof(int32_t), 1, input, NCHILADA_SWAP);
  fread_mswap(&nh->ndim, sizeof(int32_t), 1, input, NCHILADA_SWAP);
  fread_mswap(&nh->code, sizeof(int32_t), 1, input, NCHILADA_SWAP);
  if (nh->ndim != dims) {
    fprintf(stderr, "[Error] Unexpected dimensions in NChilada file %s (%d instead of %"PRId64")\n", filename, nh->ndim, dims);
    exit(1);
  }
  int64_t total_p = nh->iHighWord;
  total_p <<= 32;
  total_p += nh->nbodies;
  return total_p;
}


void read_nchilada(char *filename_stem, char *ext, struct particle **p, int64_t *num_p, int64_t offset, int64_t width, int64_t stride, int64_t expand, int64_t block) {
  char buffer[1024];
  snprintf(buffer, 1024, "%s/%s", filename_stem, ext);
  FILE *input = check_fopen(buffer, "r");
  struct nchilada_header nh;
  int64_t total_p = read_nchilada_header(input, buffer, &nh, stride);
  SCALE_NOW = nh.time;
  int64_t filetype_width = 0;
  void *rbuffer = NULL;

  filetype_width = (nh.code == float32) ? sizeof(float) : 
    (nh.code == float64) ? sizeof(double) : 
    (nh.code == int32 || nh.code == uint32) ? sizeof(int32_t) : 
    (nh.code == int64 || nh.code == uint64) ? sizeof(int64_t) : 0;
  if (!filetype_width) {
    fprintf(stderr, "[Error] Unsupported data type in file %s!\n", buffer);
    exit(1);
  }

  if (filetype_width > width && (nh.code != float64)) {
    fprintf(stderr, "[Error] Data element width in file %s is too large to fit into Rockstar's particle format.\n", buffer);
    exit(1);
  }

  if (offset==OFFSET(id)) NCHILADA_ID_BYTES = filetype_width;

  int64_t p_start, to_read;
  particle_range(total_p, block, NUM_BLOCKS, &p_start, &to_read);
  if (!(to_read>0)) { fclose(input); return; }
  p_start += 2; //Ignore minimum + maximum values at beginning of file
  check_fseeko(input, filetype_width*stride*p_start, SEEK_CUR);

  check_realloc_s(rbuffer, filetype_width*stride*NCHILADA_RBUFFER, sizeof(char));
  int64_t pp_start = *num_p;
  if (!expand) { pp_start -= to_read; }
  else {
    *num_p = *num_p + to_read;
    check_realloc_s(*p, (*num_p), sizeof(struct particle));
    memset(*p+pp_start, 0, sizeof(struct particle)*(to_read));
  }

  int64_t total_read = 0, i, n;
  while (total_read < to_read) {
    n = to_read - total_read;
    if (n > NCHILADA_RBUFFER) n = NCHILADA_RBUFFER;
    fread_mswap(rbuffer, filetype_width, stride*n, input, NCHILADA_SWAP);
    if ((stride != 3) || (width != 8)) {
      for (i=0; i<n; i++) memcpy(((char *)&(p[0][pp_start+i])) + offset,
                                 rbuffer+(i*stride*filetype_width),
				 stride*filetype_width);
    }
    else {
      for (i=0; i<n; i++) {
        float *dest = (float *)((char *)&(p[0][pp_start+i]) + offset);
        for (int64_t j=0; j<3; j++)
          dest[j] = (float)(*((double *)(rbuffer + (i*stride*filetype_width) + j*filetype_width)));
      }
    }
    pp_start += n;
    total_read += n;
  }
  free(rbuffer);
  fclose(input);
}

void convert_nchilada(int64_t p_start, struct particle *p, int64_t num_p) {
  uint32_t id;
  for (int64_t i=p_start; i<num_p; i++) {
    for (int64_t j=0; j<3; j++) {
      p[i].pos[j] = (p[i].pos[j]+0.5)*NCHILADA_LENGTH_CONVERSION;
      p[i].pos[j+3] *= NCHILADA_VELOCITY_CONVERSION *SCALE_NOW;
    }
    if (NCHILADA_ID_BYTES == 4) {
      memcpy(&id, &(p[i].id), sizeof(uint32_t));
      p[i].id = id;
    }
    p[i].mass *= NCHILADA_MASS_CONVERSION;
    p[i].energy *= ENERGY_PER_KELVIN_PER_UNIT_MASS;
    p[i].softening *= NCHILADA_LENGTH_CONVERSION;
  }
}

void load_particles_nchilada(char *filename, struct particle **p, int64_t *num_p) {
  int64_t length = strlen(filename);
  int64_t block = 0;

  if (NUM_BLOCKS > 1) {
    length--;
    while (length && filename[length] != '.') length--;
    assert(filename[length]=='.');
    filename[length] = 0;
    block = atol(filename+length+1);
  }

  int64_t gas_start = *num_p;
  //Gas:
  read_nchilada(filename, "gas/iord", p, num_p, OFFSET(id), sizeof(int64_t), 1, 1, block);
  read_nchilada(filename, "gas/pos", p, num_p, OFFSET(pos[0]), sizeof(float), 3, 0, block);
  read_nchilada(filename, "gas/vel", p, num_p, OFFSET(pos[3]), sizeof(float), 3, 0, block);
  read_nchilada(filename, "gas/mass", p, num_p, OFFSET(mass), sizeof(float), 1, 0, block);
  read_nchilada(filename, "gas/temperature", p, num_p, OFFSET(energy), sizeof(float), 1, 0, block);
  read_nchilada(filename, "gas/metals", p, num_p, OFFSET(metallicity), sizeof(float), 1, 0, block);
  read_nchilada(filename, "gas/soft", p, num_p, OFFSET(softening), sizeof(float), 1, 0, block);
  

  //Dark Matter
  int64_t dark_start = *num_p;
  read_nchilada(filename, "dark/iord", p, num_p, OFFSET(id), sizeof(int64_t), 1, 1, block);
  read_nchilada(filename, "dark/pos", p, num_p, OFFSET(pos[0]), sizeof(float), 3, 0, block);
  read_nchilada(filename, "dark/vel", p, num_p, OFFSET(pos[3]), sizeof(float), 3, 0, block);
  read_nchilada(filename, "dark/mass", p, num_p, OFFSET(mass), sizeof(float), 1, 0, block);
  read_nchilada(filename, "dark/soft", p, num_p, OFFSET(softening), sizeof(float), 1, 0, block);

  int64_t stars_start = *num_p;
  //Stars
read_nchilada(filename, "star/iord", p, num_p, OFFSET(id), sizeof(int64_t), 1, 1, block);
  read_nchilada(filename, "star/pos", p, num_p, OFFSET(pos[0]), sizeof(float), 3, 0, block);
  read_nchilada(filename, "star/vel", p, num_p, OFFSET(pos[3]), sizeof(float), 3, 0, block);
  read_nchilada(filename, "star/mass", p, num_p, OFFSET(mass), sizeof(float), 1, 0, block);
  read_nchilada(filename, "star/metals", p, num_p, OFFSET(metallicity), sizeof(float), 1, 0, block);
  read_nchilada(filename, "star/soft", p, num_p, OFFSET(softening), sizeof(float), 1, 0, block);

  if (NUM_BLOCKS > 1) filename[length] = '.';

  int64_t i;
  for (i=gas_start; i<dark_start; i++) p[0][i].type = RTYPE_GAS;
  for (; i<stars_start; i++) p[0][i].type = RTYPE_DM;
  for (; i<*num_p; i++) p[0][i].type = RTYPE_STAR;

  convert_nchilada(gas_start, *p, *num_p);
}

#undef OFFSET
