#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "check_syscalls.h"
#include "fof.h"
#include "particle.h"
#include "config_vars.h"

struct fof *fofs = NULL;
struct smallfof *smallfofs = NULL;
int64_t num_fofs=0, num_alloced_fofs = 0, 
  num_smallfofs = 0, num_alloced_smallfofs = 0;
struct particle *root_p = 0;
int64_t *particle_smallfofs = NULL;
int64_t num_particles = 0, num_alloced_particles = 0;
int64_t num_boundary_fofs = 0;

void init_particle_smallfofs(int64_t num_p, struct particle *particles) {
  int64_t i;
  if (num_p > num_alloced_particles) {
    particle_smallfofs = 
      check_realloc(particle_smallfofs, sizeof(int64_t)*num_p,
		    "Allocating particle smallfof links.");
    num_alloced_particles = num_p;
  }
  for (i=0; i<num_p; i++) particle_smallfofs[i] = -1;
  root_p = particles;
  num_particles = num_p;
  num_boundary_fofs = num_fofs = num_smallfofs = 0;
}

int64_t add_new_smallfof(void) {
  if (num_smallfofs >= num_alloced_smallfofs) {
    smallfofs = (struct smallfof *)
      check_realloc(smallfofs, sizeof(struct smallfof)*(num_smallfofs+1000),
		    "Allocating SmallFOFs.\n");
    num_alloced_smallfofs += 1000;
  }
  smallfofs[num_smallfofs].root = num_smallfofs;
  num_smallfofs++;
  return (num_smallfofs-1);
}

void _collapse_smallfof(struct smallfof *f) {
  struct smallfof *r = smallfofs + f->root, *n;
  if (smallfofs[f->root].root == f->root) return;
  while (r != (smallfofs + r->root)) r = smallfofs + r->root;
  while (f != r) {
    n = smallfofs + f->root;
    f->root = r->root;
    f = n;
  }
}

void merge_smallfofs(struct smallfof *f1, struct smallfof *f2) {
  struct smallfof *r = NULL, *d;
  int64_t f1root;
  if (f2->root == f1->root) return;
  _collapse_smallfof(f1);
  if (f2 == smallfofs + f2->root) { f2->root = f1->root; return; }
  f1root = f1->root;
  d = smallfofs + f1root;
  while (r != d) {
    r = smallfofs + f2->root;
    f2->root = f1root;
    f2 = r;
  }
}

int64_t tag_boundary_particle(struct particle *p) {
#define SMALLFOF_OF(a) particle_smallfofs[(a) - root_p]
  int64_t f = SMALLFOF_OF(p);
  if (f<0) {
    SMALLFOF_OF(p) = add_new_smallfof();
    num_boundary_fofs++;
    return num_boundary_fofs-1;
  }
  _collapse_smallfof(smallfofs + f);
  if (smallfofs[f].root >= num_smallfofs - num_boundary_fofs)
    return (smallfofs[f].root - (num_smallfofs - num_boundary_fofs));
  int64_t new_fof = add_new_smallfof();
  smallfofs[smallfofs[f].root].root = new_fof;
  num_boundary_fofs++;
  return (num_boundary_fofs-1);
#undef SMALLFOF_OF
}


void link_particle_to_fof(struct particle *p, int64_t n, struct particle **links) {
  int64_t i;
#define SMALLFOF_OF(a) particle_smallfofs[(a) - root_p]
  int64_t f = SMALLFOF_OF(p);
  if (n<2) return;
  if (f<0) {
    for (i=0; i<n; i++) {
      if (SMALLFOF_OF(links[i]) == -1) continue;
      else { f = SMALLFOF_OF(links[i]); break; }
    }
    if (f<0) { f = add_new_smallfof(); }
  }
  for (i=0; i<n; i++) {
    if (SMALLFOF_OF(links[i]) == -1) SMALLFOF_OF(links[i]) = f;
    else merge_smallfofs(smallfofs + SMALLFOF_OF(links[i]), smallfofs + f);
  }
#undef SMALLFOF_OF
}

void link_fof_to_fof(struct particle *p, int64_t n, struct particle **links) {
  int64_t i;
#define SMALLFOF_OF(a) particle_smallfofs[(a) - root_p]
  int64_t f = SMALLFOF_OF(p);
  if ((n<2) || (f<0)) return;
  for (i=0; i<n; i++) {
    if (SMALLFOF_OF(links[i]) == -1) continue;
    if (SMALLFOF_OF(links[i]) == f) continue;
    merge_smallfofs(smallfofs + SMALLFOF_OF(links[i]), smallfofs + f);
  }
#undef SMALLFOF_OF
}

void collapse_smallfofs(void) {
  int64_t i;
  for (i=0; i<num_smallfofs; i++) _collapse_smallfof(smallfofs + i);
  for (i=0; i<num_particles; i++) {
    if (particle_smallfofs[i] < 0) continue;
    particle_smallfofs[i] = smallfofs[particle_smallfofs[i]].root;
  }
}

int64_t add_new_fof(void) {
  if (num_fofs >= num_alloced_fofs) {
    fofs = (struct fof *)
      check_realloc(fofs, sizeof(struct fof)*(num_fofs+1000),
		    "Allocating FOFs.\n");
    memset(fofs + num_fofs, 0, sizeof(struct fof)*1000);
    num_alloced_fofs = num_fofs + 1000;
  }
  num_fofs++;
  return (num_fofs-1);
}

void partition_sort_particles(int64_t min, int64_t max,
		struct particle *particles, int64_t *assignments) {  
  int64_t minpivot, maxpivot, pivot, i, si, tmp;
  struct particle tmp_p;
  if (max-min < 2) return;

  maxpivot = minpivot = assignments[min];
  for (i=min+1; i<max; i++) {
    if (assignments[i] > maxpivot) maxpivot = assignments[i];
    if (assignments[i] < minpivot) minpivot = assignments[i];
  }
  if (minpivot==maxpivot) return;
  pivot = minpivot + (maxpivot-minpivot)/2;
  si = max-1;
#define SWAP(a,b) {tmp_p = particles[a]; particles[a] = particles[b]; \
    particles[b] = tmp_p; tmp = assignments[a];			      \
    assignments[a] = assignments[b];  assignments[b] = tmp;}
  
  for (i=min; i<si; i++)
    if (assignments[i] > pivot) { SWAP(i, si); si--; i--; }
  if (i==si && assignments[si]<=pivot) si++;
#undef SWAP
  partition_sort_particles(min, si, particles, assignments);
  partition_sort_particles(si, max, particles, assignments);
}

void build_fullfofs(void) {
  int64_t i, sf=-1, last_sf=-1, f=-1;
  collapse_smallfofs();
  partition_sort_particles(0, num_particles, root_p, particle_smallfofs);
  for (i=0; i<num_particles; i++) {
    if (particle_smallfofs[i] < 0) continue;
    sf = particle_smallfofs[i];
    if (sf != last_sf) {
      if (f>-1) { fofs[f].num_p = (root_p + i) - fofs[f].particles; }
      if ((f < 0) || (fofs[f].num_p >= MIN_HALO_PARTICLES) ||
	  last_sf >= num_smallfofs-num_boundary_fofs) f=add_new_fof();
      fofs[f].particles = root_p + i;
      last_sf = sf;
    }
  }
  if (f>-1) { 
    fofs[f].num_p = (root_p + i) - fofs[f].particles;
    if (fofs[f].num_p < MIN_HALO_PARTICLES && 
	(sf < num_smallfofs-num_boundary_fofs)) num_fofs--;
  }
  num_smallfofs = 0;
}

struct fof *return_fullfofs(int64_t *num_f, int64_t *num_bf) {
  struct fof *f_all_fofs;
  num_smallfofs = num_alloced_smallfofs = 0;
  smallfofs = check_realloc(smallfofs, 0, "Freeing SmallFOFs.");
  particle_smallfofs = check_realloc(particle_smallfofs, 0,
				     "Freeing particle smallfofs.");
  num_alloced_particles = 0;
  f_all_fofs = fofs;
  *num_f = num_fofs;
  *num_bf = num_boundary_fofs;
  fofs = NULL;
  num_boundary_fofs = num_alloced_fofs = num_fofs = 0;
  return f_all_fofs;
}


void copy_fullfofs(struct fof **base, int64_t *num_f, int64_t *num_alloced_f)
{
  if ((*num_f)+num_fofs > (*num_alloced_f)) {
    *base = check_realloc(*base, sizeof(struct fof)*((*num_f)+num_fofs+1000),
			  "Allocating copy space for FOFs.");
    *num_alloced_f = (*num_f)+num_fofs+1000;
  }
  memcpy((*base)+(*num_f), fofs, sizeof(struct fof)*num_fofs);
  *num_f = (*num_f) + num_fofs;

  //Make sure to delete last fof, if necessary.
  if (num_fofs < num_alloced_fofs) num_fofs++;
  memset(fofs, 0, sizeof(struct fof)*num_fofs);
  num_boundary_fofs = num_fofs = 0;
}
