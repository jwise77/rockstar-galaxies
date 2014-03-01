#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <string.h>
#include <assert.h>
#include "../config_vars.h"
#include "../config.h"
#include "../check_syscalls.h"
#include "../rockstar.h"
#include "../io/meta_io.h"
#include "../groupies.h"

#define FAST3TREE_TYPE struct halo
#define FAST3TREE_PREFIX SUBSTATS
#include "../fast3tree.c"
#define RADIUS r
#define parent desc
#define GROUP_LIST halos
#include "../parents.c"

int sort_halos_by_id(const void *a, const void *b) {
  const struct halo *c = a;
  const struct halo *d = b;
  if (c->id < d->id) return -1;
  if (c->id > d->id) return 1;
  return 0;
}

double dot(float *x1, float *x2) {
  int64_t i;
  double s=0;
  for (i=0; i<3; i++) s+=x1[i]*x2[i];
  return s;
}

void cross(float *x1, float *x2, float *x3) {
  x3[0] = x3[1] = x3[2] = 0;
#define cross(a,x,y,s) x3[a] s (x1[x]*x2[y])
  cross(0,1,2,+=);
  cross(0,2,1,-=);
  cross(1,2,0,+=);
  cross(1,0,2,-=);
  cross(2,0,1,+=);
  cross(2,1,0,-=);
#undef cross
}

float calc_angle(float *pos1, float *pos2) {
  double s1=0, s2=0, d=0, norm;
  s1 = dot(pos1, pos1);
  s2 = dot(pos2, pos2);
  d = dot(pos1, pos2);
  norm = sqrt(s1*s2);
  if (!norm) return 0;
  return acos(d/norm);
}



void decompose_v(float *pos, float *axis, float *vel, float *vr, float *vtheta, float *vphi) {
  float theta[3], phi[3];
  float tnorm, pnorm;
  float r = sqrt(dot(pos, pos));
  *vr = *vtheta = *vphi = 0;
  if (!r) return;
  *vr = dot(vel, pos) / r;

  cross(axis, pos, theta);
  tnorm = sqrt(dot(theta, theta));
  if (!tnorm) return;
  *vtheta = dot(vel, theta) / tnorm;

  cross(pos, theta, phi);
  pnorm = sqrt(dot(phi, phi));
  *vphi = dot(vel, phi) / pnorm;
}

void _calc_subhalo_stats(struct halo *h, FILE *output)
{
  struct halo *parent = halos + h->parent;
  float pos[6], r=0, theta, vr, vtheta, vphi;
  int64_t k;

  while (parent->parent >= 0) parent = halos + parent->parent;
  if (!parent->m || !parent->r || !parent->vmax || !parent->rs) return;

  for (k=0; k<6; k++) pos[k] = h->pos[k] - parent->pos[k];
  for (k=0; k<3; k++) r += pos[k]*pos[k];
  r = sqrt(r);
  theta = calc_angle(pos, parent->J);
  decompose_v(pos, parent->J, pos+3, &vr, &vtheta, &vphi);
  fprintf(output, "%"PRId64" %"PRId64" %.3e %.3e %.3f %.3f %.3f %.3f %.1f "
	  "%.3f %.3f %.3f %.3f %.3f\n", h->id, parent->id, h->m/parent->m,
	  parent->m, h->vmax/parent->vmax, parent->vmax, h->r/parent->r, 
	  parent->r/1e3, parent->r / parent->rs, r*1e3/parent->r, theta,
	  vr/parent->vmax, vtheta/parent->vmax, vphi/parent->vmax);
}

void calc_subhalo_stats(int64_t snap)
{
  int64_t i;
  FILE *input, *output;
  struct binary_output_header hdr;
  char buffer[1024];
  num_halos = 0;
  for (i=0; i<NUM_WRITERS; i++) {
    get_output_filename(buffer, 1024, snap, i, "bin");
    input = check_fopen(buffer, "rb");
    check_fread(&hdr, sizeof(struct binary_output_header), 1, input);
    halos = check_realloc(halos, sizeof(struct halo)*(num_halos+hdr.num_halos),
			  "Allocating room for halos.");
    check_fread(halos + num_halos, sizeof(struct halo), hdr.num_halos, input);
    num_halos += hdr.num_halos;
    fclose(input);
  }

  find_parents(num_halos);
  qsort(halos, num_halos, sizeof(struct halo), sort_halos_by_id);
  get_output_filename(buffer, 1024, snap, 0, "subs");
  output = check_fopen(buffer, "w");
  fprintf(output, "#ID UPID Mass/Pmass Pmass Vmax/Pvmax Pvmax Rvir/Prvir Prvir Pc R/Rvir Theta Vr/Pvmax Vtheta/Pvmax Vphi/Pvmax\n");
  for (i=0; i<num_halos; i++) {
    if (halos[i].parent < 0) continue;
    _calc_subhalo_stats(halos+i, output);
  }
  fclose(output);
}

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

  calc_subhalo_stats(snap);
  return 0;
}


