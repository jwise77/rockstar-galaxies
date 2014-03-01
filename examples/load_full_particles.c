/* Example code for reading in full particle dumps. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <math.h>
#include <assert.h>
#include <sys/stat.h>
#include "../check_syscalls.h"
#include "../io/stringparse.h"
#include "../halo.h"
#include "../particle.h"
#include "../config_vars.h"
#include "load_full_particles.h"


#ifndef TEST_LOADFP
extern struct halo *h; //(Must call list of halos "h"!)
#else
struct halo *h = NULL;
struct full_particle *p = NULL;
float bounds[6] = {0};
int64_t num_h=0, num_p=0;

int main(int argc, char **argv) {
  int64_t i;

  if (argc < 2) {
    printf("Usage: %s file1.particles ...\n", argv[0]);
    exit(1);
  }

  for (i=1; i<argc; i++) {
    num_h = num_p = 0;
    load_full_particles(argv[i], &h, &num_h, &p, &num_p, bounds);
    printf("Scale factor: %f\n", SCALE_NOW);
    printf("Bounds: (%f, %f, %f) - (%f, %f, %f)\n", bounds[0], bounds[1], 
	   bounds[2], bounds[3], bounds[4], bounds[5]);
    printf("Om = %f; Ol = %f; h0 = %f\n", Om, Ol, h0);
    printf("Particle mass: %e Msun/h\n", PARTICLE_MASS);
    printf("Num groups: %"PRId64"\n", num_h);
    printf("Num particles: %"PRId64"\n", num_p);
    //Do stuff with particles here...
    //See io/bgc2.h for a description of the group and particle structure format
  }
  return 0;
}
#endif /* TEST_LOADFP */

#define FAST3TREE_TYPE struct halo
#define IGNORE_INVALID_IDS
#define NO_PERIODIC_BOUNDARIES
#include "../fast3tree.c"
#include "../parents.c"

int int_id_compare(const void *a, const void *b) {
  const struct halo *c = a;
  const struct halo *d = b;
  if (c->flags < d->flags) return -1;
  if (c->flags > d->flags) return 1;
  return 0;
}


#define HALO_FIELDS 20
#define PARTICLE_FIELDS 13

void load_full_particles(char *filename, struct halo **h, int64_t *num_h, 
			 struct full_particle **p, int64_t *num_p, float *bnds) {
  FILE *input;
  int64_t i;
  char buffer[1024];
  struct halo the_h = {0};
  struct full_particle the_p;
  SHORT_PARSETYPE;
  enum short_parsetype stypes_h[HALO_FIELDS] = 
    { D64, D64, D64, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F};
  enum parsetype types_h[HALO_FIELDS];
  enum short_parsetype stypes_p[PARTICLE_FIELDS] = 
    { F, F, F, F, F, F, F, F, D64, D, D64, D64, D64};
  enum parsetype types_p[PARTICLE_FIELDS];

  void *data_h[HALO_FIELDS] = {&(the_h.id), &(the_h.flags), &(the_h.num_p), &(the_h.m),
		      &(the_h.mgrav), &(the_h.r), &(the_h.vmax), &(the_h.rvmax),
		      &(the_h.vrms), the_h.pos, the_h.pos+1, the_h.pos+2,
		      the_h.pos+3, the_h.pos+4, the_h.pos+5, the_h.J,
		      the_h.J+1, the_h.J+2, &(the_h.energy), &(the_h.spin)};
  void *data_p[PARTICLE_FIELDS] = {the_p.pos, the_p.pos+1, the_p.pos+2, 
		      the_p.pos+3, the_p.pos+4, the_p.pos+5, &(the_p.mass),
		      &(the_p.energy), &(the_p.id), &(the_p.type), &(the_p.a_hid),
		      &(the_p.hid), &(the_p.ehid)};

  input = check_fopen(filename, "r");
  for (i=0; i<HALO_FIELDS; i++) types_h[i] = stypes_h[i];
  for (i=0; i<PARTICLE_FIELDS; i++) types_p[i] = stypes_p[i];

  //Header info
  while (fgets(buffer, 1024, input)) {
    if (!strcmp("#Halo table begins here:\n", buffer)) break;
    if (!strncmp(buffer, "#a = ", 5)) SCALE_NOW = atof(buffer+5);
    if (!strncmp(buffer, "#Om = ", 6))
      sscanf(buffer, "#Om = %lf; Ol = %lf; h = %lf", &Om, &Ol, &h0);
    if (!strncmp(buffer, "#Particle mass: ", 16))
      PARTICLE_MASS = atof(buffer+16);
    if (!strncmp(buffer, "#Bounds: ", 9))
      sscanf(buffer, "#Bounds: (%f, %f, %f) - (%f, %f, %f)",
	     bnds, bnds+1, bnds+2, bnds+3, bnds+4, bnds+5);
    if (!strncmp(buffer, "#Box size: ", 11))
      BOX_SIZE = atof(buffer+11);
  }

  //Halos
  while (fgets(buffer, 1024, input)) {
    if (!strcmp("#Particle table begins here:\n", buffer)) break;
    i = stringparse(buffer+1, data_h, (enum parsetype *)types_h, HALO_FIELDS);
    if (i < HALO_FIELDS) continue;
    if (!((*num_h)%1000)) {
      *h = check_realloc(*h, sizeof(struct halo)*((*num_h)+1000),
			 "Allocating halos");
    }
    (*h)[*num_h] = the_h;
    *num_h = *num_h + 1;
  }

  find_parents(*num_h);
  qsort(*h, *num_h, sizeof(struct halo), int_id_compare);

  //Particles
  while (fgets(buffer, 1024, input)) {
    i = stringparse(buffer, data_p, (enum parsetype *)types_p, PARTICLE_FIELDS);
    if (i < PARTICLE_FIELDS) continue;
    if (!((*num_p)%1000)) {
      *p = check_realloc(*p, sizeof(struct full_particle)*((*num_p)+1000),
			 "Allocating particles");
    }
    (*p)[*num_p] = the_p;
    *num_p = *num_p + 1;
  }
  fclose(input);
}


