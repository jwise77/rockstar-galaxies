#include <stdio.h>
#include <stdlib.h>
#include <sys/mman.h>
#include "io/meta_io.h"
#include "halo.h"
#include "particle.h"
#include "config_vars.h"
#include "check_syscalls.h"
#include "universe_time.h"
#include "inthash.h"
#include "rockstar.h"
#include "bounds.h"
#include "fun_times.h"

#define DEBUG_FUN_TIMES 0

#define FAST3TREE_TYPE struct previous_halo
#define FAST3TREE_PREFIX FUN_TIMES
#include "fast3tree.c"

struct previous_halo *ph = NULL;
int64_t prev_snap = -1;
int64_t num_prev_halos = 0;
struct halo *prev_halo_buffer = NULL;
struct fast3tree *phtree = NULL;
struct fast3tree_results *phtree_res = NULL;
void **prev_files = NULL;
int64_t *prev_file_lengths = NULL;
int64_t *prev_chunks = NULL;
int64_t num_prev_files = 0;
int64_t max_num_p = 0;
struct prev_bounds *p_bounds = NULL;


static inline void add_to_previous_halos(struct halo *h, struct binary_output_header *bh, void *file) {
  if (!(num_prev_halos%1000))
    check_realloc_s(ph, sizeof(struct previous_halo), (num_prev_halos+1000));
  memcpy(ph[num_prev_halos].pos, h->pos, sizeof(float)*6);
  ph[num_prev_halos].m = h->m;
  ph[num_prev_halos].r = h->r;
  ph[num_prev_halos].id = h->id;
  ph[num_prev_halos].num_p = h->num_p;
  ph[num_prev_halos].chunk = bh->chunk;
  ph[num_prev_halos].particles = (void *)(file + sizeof(struct binary_output_header) +
     sizeof(struct halo)*bh->num_halos + h->p_start*sizeof(int64_t));
  num_prev_halos++;
}

void *add_new_prev_file(char *filename, int64_t chunk) {
  for (int64_t i=0; i<num_prev_files; i++)
    if (prev_chunks[i]==chunk) return NULL;
  if (!(num_prev_files % 10)) {
    check_realloc_s(prev_files, sizeof(void *),(num_prev_files+10));
    check_realloc_s(prev_file_lengths, sizeof(int64_t),(num_prev_files+10));
    check_realloc_s(prev_chunks, sizeof(int64_t),(num_prev_files+10));
  }
  prev_files[num_prev_files] = 
    check_mmap_file(filename, 'r', prev_file_lengths + num_prev_files);
  prev_chunks[num_prev_files] = chunk;
  num_prev_files++;
  return prev_files[num_prev_files-1];
}

void clear_prev_files(void) {
  for (int64_t i=0; i<num_prev_files; i++)
    munmap(prev_files[i], prev_file_lengths[i]);
  num_prev_files = 0;
  num_prev_halos = 0;
}

void load_prev_binary_halos(int64_t snap, int64_t chunk, float *bounds, int64_t our_chunk) {
  void *input;
  char buffer[1024];
  struct binary_output_header bh;
  int64_t offset, remaining = 0, to_read, i,j;
  double v_to_dx;

  max_num_p = 0;
  get_output_filename(buffer, 1024, snap, chunk, "bin");
  input = add_new_prev_file(buffer, chunk);
  if (!input) return; //Already loaded file
  memcpy(&bh, input, sizeof(struct binary_output_header));
  assert(bh.magic == ROCKSTAR_MAGIC);
  assert(bh.num_halos >= 0);
  assert(bh.num_particles >= 0);

  //Conversion in Comoving Mpc/h / (km/s)
  //Note that the time units are in 1/H = 1/(h*100 km/s/Mpc)
  v_to_dx = 0.01*(scale_to_time(SCALE_NOW) - scale_to_time(bh.scale)) / 
    (0.5*(SCALE_NOW + bh.scale));

  remaining = bh.num_halos;
  offset = sizeof(struct binary_output_header);
  while (remaining > 0) {
    to_read = PREV_HALO_BUFFER_SIZE;
    if (to_read > remaining) to_read = remaining;
    memcpy(prev_halo_buffer, input+offset, sizeof(struct halo)*to_read);
    remaining -= to_read;
    offset += sizeof(struct halo)*to_read;
    for (i=0; i<to_read; i++) {
      for (j=0; j<3; j++) 
	prev_halo_buffer[i].pos[j] += v_to_dx*prev_halo_buffer[i].pos[j+3];
      add_to_previous_halos(prev_halo_buffer + i, &bh, input);
    }
  }
}

void load_previous_halos(int64_t snap, int64_t chunk, float *bounds) {
  int64_t rchunk;
  struct binary_output_header bh;
  float overlap_region[6];

  if (!phtree) {
    phtree = fast3tree_init(num_prev_halos, ph);
    phtree_res = fast3tree_results_init();
  }
  if (snap == STARTING_SNAP || LIGHTCONE || !PARALLEL_IO) return;
  check_realloc_s(prev_halo_buffer,sizeof(struct halo),PREV_HALO_BUFFER_SIZE);
  if (prev_snap != (snap-1)) {
    if (!p_bounds)
      check_realloc_s(p_bounds, sizeof(struct prev_bounds), NUM_WRITERS);
    prev_snap = snap-1;
    for (rchunk = 0; rchunk < NUM_WRITERS; rchunk++) {
      load_binary_header(snap-1, rchunk, &bh);
      memcpy(p_bounds[rchunk].bounds, bh.bounds, sizeof(float)*6);
    }
  }

  for (rchunk = 0; rchunk < NUM_WRITERS; rchunk++) {
    if (!bounds || bounds_overlap(p_bounds[rchunk].bounds, bounds,
				  overlap_region, OVERLAP_LENGTH))
      load_prev_binary_halos(snap-1, rchunk, bounds, chunk);
  }
  check_realloc_s(prev_halo_buffer, 0, 0);
  fast3tree_rebuild(phtree, num_prev_halos, ph);
}

int64_t prev_halo_acceptable(struct halo *h, struct previous_halo *tph) {
  int64_t j;
  float ds, dx=0, dv=0;
  if (tph->num_p < h->num_p*0.5) return 0;
  for (j=0; j<6; j++) {
    ds = h->pos[j]-tph->pos[j];
    if (j < 3) dx += ds*ds;
    else dv += ds*ds;
  }
  if (dv > h->vrms*h->vrms) return 0;
  return 1;
}

int sort_previous_halos(const void *a, const void *b) {
  struct previous_halo *c = *((struct previous_halo **)a);
  struct previous_halo *d = *((struct previous_halo **)b);
  if (c->particles < d->particles) return -1;
  if (c->particles > d->particles) return 1;
  return 0;
}

int sort_core_particles(const void *a, const void *b) {
  const struct particle *c = a;
  const struct particle *d = b;
  if (c->pos[0] < d->pos[0]) return -1;
  if (c->pos[0] > d->pos[0]) return 1;
  return 0;
}

void convert_and_sort_core_particles(struct halo *h, struct particle *hp, float max_r, int64_t *n_core) {
  int64_t i, j, np=h->num_p;
  float ds, dx;
  for (i=0; i<np; i++) {
    for (j=0,dx=0; j<3; j++) { ds = h->pos[j]-hp[i].pos[j]; dx+= ds*ds; }
    hp[i].pos[0] = ds;
    if (max_r && (ds > max_r * max_r)) {
      struct particle tmp = hp[i];
      hp[i] = hp[np-1];
      hp[np-1] = tmp;
      np--;
      i--;
    }
  }
  qsort(hp, np, sizeof(struct particle), sort_core_particles);
  if (n_core) *n_core = np;
  for (i=0; i<h->num_p; i++)
    hp[i].pos[0] = p[hp[i].id].pos[0];
}

void reassign_particles_to_parent(struct halo *h, struct particle *hp, int64_t *particle_halos, struct halo *parent_h) {
  int64_t i, prev_id = extra_info[parent_h-halos].ph;
  struct previous_halo *tph=NULL;
  struct inthash *ih = NULL;
  if (prev_id < 0) return;
  fast3tree_find_sphere(phtree, phtree_res, parent_h->pos, parent_h->r);
  if (!phtree_res->num_points) return;
  for (i=0; i<phtree_res->num_points; i++) {
    if (phtree_res->points[i]->id == prev_id) {
      tph = phtree_res->points[i];
      break;
    }    
  }
  if (!tph) return;

  ih = new_inthash();
  mlock(tph->particles, tph->num_p*sizeof(int64_t));
  for (i=0; i<tph->num_p; i++)
    ih_setval(ih, tph->particles[i], (void *)1);
  munlock(tph->particles, tph->num_p*sizeof(int64_t));

  for (i=0; i<h->num_p; i++)
    if (ih_getval(ih,  p[hp[i].id].id))
      particle_halos[i] = parent_h - halos;

  free_inthash(ih);
  extra_info[h-halos].ph = -1;
}

float find_previous_mass(struct halo *h, struct particle *hp, int64_t *best_num_p, float max_r) {
  int64_t i, j, max_particles, best_particles = 0, cur_part, remaining;
  struct previous_halo *tph, *best_ph=NULL;
  struct inthash *ih = NULL;

  *best_num_p = 0;
  if (h->num_p < 100 || !num_prev_halos) return 0;
  fast3tree_find_sphere(phtree, phtree_res, h->pos, h->r);
  if (!phtree_res->num_points) return 0;
  for (i=0; i<phtree_res->num_points; i++) {
    if (!prev_halo_acceptable(h, phtree_res->points[i])) {
      phtree_res->num_points--;
      phtree_res->points[i] = phtree_res->points[phtree_res->num_points];
      i--;
    }
  }
  if (!phtree_res->num_points) return 0;

  qsort(phtree_res->points, phtree_res->num_points, sizeof(struct previous_halo *),
	sort_previous_halos);

  //convert_and_sort_core_particles(h, hp, 100.0*max_r, &max_particles);
  max_particles = h->num_p;

  ih = new_inthash();
  for (i=0; i<max_particles; i++)
    ih_setval(ih, p[hp[i].id].id, (void *)1);

  for (i=0; i<phtree_res->num_points; i++) {
    tph = phtree_res->points[i];
    remaining = tph->num_p;
    cur_part = 0;
    /*if (remaining*0.1 > MAX_CORE_PARTICLES) remaining *= 0.1;
    else if (remaining > MAX_CORE_PARTICLES) remaining = MAX_CORE_PARTICLES;
    */

    mlock(tph->particles, tph->num_p*sizeof(int64_t));
    for (j=0; j<remaining; j++)
      if (ih_getval(ih, tph->particles[j])) cur_part++;
    munlock(tph->particles, tph->num_p*sizeof(int64_t));

    if (cur_part > best_particles) {
      best_particles = cur_part;
      best_ph = tph;
    }
  }

  free_inthash(ih);
  if (best_ph) extra_info[h-halos].ph = best_ph->id;
#if DEBUG_FUN_TIMES
  if (best_ph && best_ph->m > 1e13 && best_particles > max_particles*0.1) {
    fprintf(stderr, "Hnow: %f %f %f (#%"PRId64"; %"PRId64"; %e); Phalo: %f %f %f (%e); %"PRId64"\n",
	    h->pos[0], h->pos[1], h->pos[2], (int64_t)(h-halos), h->num_p, h->num_p*PARTICLE_MASS,
	    best_ph->pos[0], best_ph->pos[1], best_ph->pos[2], best_ph->m,
	    best_particles);
  }
#endif /* DEBUG_FUN_TIMES */
  *best_num_p = best_particles;
  if (best_ph && (best_particles > max_particles*0.1)) return best_ph->m;
  return 0;
}
