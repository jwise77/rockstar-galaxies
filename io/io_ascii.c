#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <string.h>
#include "../check_syscalls.h"
#include "stringparse.h"
#include "../particle.h"
#include "../config_vars.h"

void gzip_file(char *filename) {
  char buffer[1024];
  snprintf(buffer, 1024, "gzip -f \"%s\"", filename);
  if (system(buffer) != 0)
    fprintf(stderr, "[Warning] Gzip of file %s failed.\n", filename);
}


void load_particles(char *filename, struct particle **p, int64_t *num_p) {
  FILE *input;
  char buffer[1024];
  int64_t n;
  struct particle d = {0};
  SHORT_PARSETYPE;
#define NUM_INPUTS 10
  enum short_parsetype stypes[NUM_INPUTS] = 
    { F, F, F, F, F, F, F, F, D64, D};
  enum parsetype types[NUM_INPUTS];
  void *data[NUM_INPUTS] = {&(d.pos[0]), &(d.pos[1]), &(d.pos[2]), &(d.pos[3]), &(d.pos[4]), &(d.pos[5]), &(d.mass), &(d.energy), &(d.id), &(d.type)};

  for (n=0; n<NUM_INPUTS; n++) types[n] = (enum parsetype)stypes[n];
  
  input = check_fopen(filename, "r");
  while (fgets(buffer, 1024, input)) {
    if (buffer[0] == '#') {
      if (!strncmp(buffer, "#a = ", 5)) SCALE_NOW = atof(buffer+5);
      continue;
    }
    n = stringparse(buffer, data, (enum parsetype *)types, NUM_INPUTS);

    if (n < NUM_INPUTS) continue;
    if (((*num_p)%1000)==0) {
      *p = check_realloc(*p, ((*num_p)+1000)*sizeof(struct particle),
			 "Adding new particles.");
    }
    //d.type = (d.energy>0) ? RTYPE_GAS : RTYPE_STAR;
    (*p)[(*num_p)] = d;
    (*num_p)++;
  }
  fclose(input);
}
