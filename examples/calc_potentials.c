#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include "../halo.h"
#include "../particle.h"
#include "../potential.h"
#include "../check_syscalls.h"
#include "load_full_particles.h"

#ifndef CALC_POTENTIALS
#error Need to define CALC_POTENTIALS when compiling this file!
#endif /* !def CALC_POTENTIALS */

struct halo *h = NULL;
int64_t *parents = NULL;
struct full_particle *p = NULL;
float bounds[6] = {0};
int64_t num_h=0, num_p=0;
struct potential *po = NULL;
int64_t max_num_po = 0;

void calc_potentials(void) {
  int64_t i,j=0,count=0,last_id=-1;
  for (i=0; i<num_p; i++) {
    if (p[i].hid != last_id) {
      if (last_id >= 0) h[last_id].num_p = i-h[last_id].p_start;
      last_id = p[i].hid;
      h[last_id].p_start = i;      
    }
  }
  if (last_id >= 0) h[last_id].num_p = i-h[last_id].p_start;

  for (i=0; i<num_h; i++) {
    struct halo *the_h = h+i;
    if (the_h->id < 0) continue;
    if (max_num_po < the_h->num_p) {
      max_num_po = the_h->num_p;
      po = check_realloc(po, sizeof(struct potential)*max_num_po, 
			 "Allocating potentials");
    }

    for (j=0; j<the_h->num_p; j++) {
      po[j].id = p[the_h->p_start+j].id;
      memcpy(po[j].pos, p[the_h->p_start+j].pos, sizeof(float)*6);
      po[j].pe = po[j].ke = 0;
    }
    compute_kinetic_energy(po, the_h->num_p, the_h->pos+3, the_h->pos);
    compute_potential(po, the_h->num_p);
    count=0;
    for (j=0; j<the_h->num_p; j++) if (po[j].pe >= po[j].ke) count++;
    printf("%"PRId64" %"PRId64"\n", the_h->id, count);
    for (j=0; j<the_h->num_p; j++) 
      if (po[j].pe >= po[j].ke) printf("%"PRId64"\n", po[j].id);
  }
}

int main(int argc, char **argv) {
  int64_t i;
  if (argc < 2) {
    printf("Usage: %s file1.particles ...\n", argv[0]);
    exit(1);
  }
  
  for (i=1; i<argc; i++) {
    num_h = num_p = 0;
    load_full_particles(argv[i], &h, &num_h, &p, &num_p, bounds);
    calc_potentials();
  }
}

