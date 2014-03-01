#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <strings.h>
#include <inttypes.h>
#include "rockstar.h"
#include "halo.h"
#include "fof.h"
#include "particle.h"
#include "groupies.h"
#include "subhalo_metric.h"
#include "check_syscalls.h"
#include "config_vars.h"
#include "universal_constants.h"
#include "potential.h"
#include "nfw.h"
#include "distance.h"
#include "fun_times.h"
#include "jacobi.h"
#include "hubble.h"

#define FAST3TREE_DIM 6
#define POINTS_PER_LEAF 40
#define FAST3TREE_PREFIX GROUPIES
#define FAST3TREE_TYPE struct particle
#include "fast3tree.c"

struct particle *copies = NULL; //For storing phase-space FOFs
int64_t *particle_halos = NULL;
float *particle_r = NULL;
struct potential *po = NULL;
int64_t num_alloc_pc = 0, num_copies = 0;

struct fof *subfofs = NULL;
int64_t num_subfofs = 0, num_alloced_subfofs = 0;

int64_t num_halos = 0;
struct halo *halos = NULL;
struct extra_halo_info *extra_info = NULL;

struct fast3tree_results *res = NULL;
struct fast3tree *phasetree = NULL;

int64_t num_alloc_gh = 0, num_growing_halos = 0;
struct halo **growing_halos = NULL;

int64_t *halo_ids = NULL;
int64_t num_alloced_halo_ids = 0;

double particle_thresh_dens[5] = {0}, particle_rvir_dens = 0,
  particle_rvir_dens_z0 = 0;
int64_t min_dens_index = 0;
double dynamical_time = 0;

float dist(float x1[3], float x2[3]) {
  double ds, dx=0;
  for (int64_t i=0; i<3; i++) { ds = x1[i]-x2[i]; dx+=ds*ds; }
  return sqrt(dx);
}

double vir_density(double a) {
  double x = (Om/pow(a,3))/pow(hubble_scaling(1.0/a-1.0),2.0) - 1.0;
  return ((18*M_PI*M_PI + 82.0*x - 39*x*x)/(1.0+x));
}

float _calc_mass_definition(char **md) {
  int64_t length = strlen(*md);
  char last_char = (length) ? md[0][length-1] : 0;
  float matter_fraction = (Om/pow(SCALE_NOW,3))/pow(hubble_scaling(1.0/SCALE_NOW-1.0),2.0);
  float cons = Om * CRITICAL_DENSITY / PARTICLE_MASS; // background density
  char *mass = *md;
  float thresh_dens;
  if (mass[0] == 'm' || mass[0] == 'M') mass++;

  if (last_char == 'b' || last_char == 'B')
    thresh_dens = atof(mass) * cons;
  else if (last_char == 'c' || last_char == 'C')
    thresh_dens = atof(mass) * cons / matter_fraction;
  else {
    if (strcasecmp(*md, "vir")) *md = "vir";
    thresh_dens = vir_density(SCALE_NOW) * cons;
  }
  particle_rvir_dens_z0 = vir_density(1.0) * cons;
  return thresh_dens;
}

void calc_mass_definition(void) {
  char *vir = "vir";
  int64_t i;
  particle_thresh_dens[0] = _calc_mass_definition(&MASS_DEFINITION);
  particle_thresh_dens[1] = _calc_mass_definition(&MASS_DEFINITION2);
  particle_thresh_dens[2] = _calc_mass_definition(&MASS_DEFINITION3);
  particle_thresh_dens[3] = _calc_mass_definition(&MASS_DEFINITION4);
  particle_thresh_dens[4] = _calc_mass_definition(&MASS_DEFINITION5);
  particle_rvir_dens = _calc_mass_definition(&vir);
  dynamical_time = 1.0/sqrt((4.0*M_PI*Gc/3.0)*particle_rvir_dens*PARTICLE_MASS);
  min_dens_index = 0;
  for (i=1; i<5; i++) 
    if (particle_thresh_dens[i] < particle_thresh_dens[min_dens_index])
      min_dens_index = i;
}

void lightcone_set_scale(float *pos) {
  int64_t i;
  float ds=0, dx=0, z;
  for (i=0; i<3; i++) { ds = pos[i]-LIGHTCONE_ORIGIN[i]; dx+=ds*ds; }
  z = comoving_distance_h_to_redshift(sqrt(dx));
  SCALE_NOW = scale_factor(z);
  calc_mass_definition();
}


void add_new_halo(void) {
  int64_t i;
  if ((num_halos % 1000)==0) {
    halos = check_realloc(halos, sizeof(struct halo)*(num_halos+1000),
			  "Allocating room for halos.");
    extra_info = check_realloc(extra_info, sizeof(struct extra_halo_info)*(num_halos+1000),
			  "Allocating room for extra halo info.");
    memset(halos+num_halos, 0, sizeof(struct halo)*1000);
    for (i=num_halos; i<num_halos+1000; i++) {
      extra_info[i].child = extra_info[i].next_cochild = 
	extra_info[i].ph = extra_info[i].prev_cochild = extra_info[i].sub_of = -1;
      extra_info[i].max_metric = 0;
      halos[i].flags |= GROWING_FLAG;
    }
  }
  num_halos++;
}

void add_more_halo_ids(void) {
  num_alloced_halo_ids+=1000;
  halo_ids = check_realloc(halo_ids, sizeof(int64_t)*num_alloced_halo_ids,
			   "Allocating room for halo ids.");
}

void add_more_growing_halos(void) {
  num_alloc_gh += 1000;
  growing_halos = check_realloc(growing_halos, sizeof(struct halo *)
				*num_alloc_gh,
				"Allocating room for growing halos.");
}

int dist_compare(const void *a, const void *b) {
  float c = ((struct potential *)a)->r2;
  float d = ((struct potential *)b)->r2;
  if (c>d) return 1;
  if (c<d) return -1;
  return 0;
}

void _reset_potentials(struct halo *base_h, struct halo *h, float *cen, int64_t p_start, int64_t level, int64_t potential_only) {
  int64_t j, k;
  float dx, r2;
  memset(po + p_start, 0, sizeof(struct potential)*h->num_p);
  for (j=0; j<h->num_p; j++) {
    r2 = 0;
    for (k=0; k<3; k++) { dx=copies[h->p_start+j].pos[k] - cen[k]; r2+=dx*dx; }
    po[p_start + j].r2 = r2;
    memcpy(po[p_start+j].pos, copies[h->p_start+j].pos, sizeof(float)*6);
    po[p_start+j].mass = copies[h->p_start+j].mass;
    po[p_start+j].energy = copies[h->p_start+j].energy;
    po[p_start+j].type = copies[h->p_start+j].type;
    if (potential_only) po[p_start + j].ke = -1;
    if (h==base_h) po[p_start + j].flags = 1;
    if (!potential_only && (h->num_p < base_h->num_p*0.03))
      po[p_start+j].flags = 2;
  }
}

int64_t calc_particle_radii(struct halo *base_h, struct halo *h, float *cen, int64_t p_start, int64_t level, int64_t potential_only) {
  int64_t j, total_p = p_start, child, first_child, parent;

  //Break accidental graph loops
  if (level >= num_alloced_halo_ids) add_more_halo_ids();
  halo_ids[level] = h-halos;
  for (j=0; j<level; j++) if (halo_ids[j] == halo_ids[level]) return p_start;

  _reset_potentials(base_h, h, cen, p_start, level, potential_only);

  first_child = child = extra_info[h-halos].child;
  total_p += h->num_p;
  while (child > -1) {
    total_p = calc_particle_radii(base_h, halos + child,
                                  cen, total_p, level+1, potential_only);
    child = extra_info[child].next_cochild;
    assert(child != first_child);
  }

  parent = extra_info[h-halos].sub_of;
  if ((h == base_h) && (parent > -1) &&
    (halos[parent].num_child_particles*INCLUDE_HOST_POTENTIAL_RATIO < h->num_child_particles)){
    total_p = calc_particle_radii(base_h, halos + parent,
                                  cen, total_p, level+1, 1);
  }
  return total_p;
}

#include "properties.c"

void _find_subfofs_at_r(struct fof *f, float target_r) {
  int64_t i;
  init_particle_smallfofs(f->num_p, f->particles);
  for (i=0; i<f->num_p; i++) {
    //if (f->particles[i].type != RTYPE_DM) continue;
    fast3tree_find_sphere(phasetree, res, f->particles[i].pos, target_r); //*pow(f->particles[i].mass/PARTICLE_MASS, 1.0/6.0));
    link_particle_to_fof(f->particles + i, res->num_points, res->points);
  }
  build_fullfofs();
}

int64_t separate_dm(struct fof *f) {
  struct particle tmp;
  int64_t num_dm = f->num_p;
  for (int64_t i=0; i<num_dm; i++) {
    if (f->particles[i].type!=RTYPE_DM) {
      num_dm--;
      tmp = f->particles[i];
      f->particles[i] = f->particles[num_dm];
      f->particles[num_dm] = tmp;
      i--;
    }
  }
  return num_dm;
}


#define MAX_PARTICLES_TO_SAMPLE 10000
void _find_subfofs_better2(struct fof *f,  float thresh) {
  int64_t i, j, num_test = MAX_PARTICLES_TO_SAMPLE;
  float target_r = 0;
  norm_sd(f, thresh);
  int64_t num_dm = separate_dm(f);
  if (!num_dm) num_dm = f->num_p;
  fast3tree_rebuild(phasetree, num_dm, f->particles);
  if (EXACT_LL_CALC) num_test = num_dm;
  if (num_test > num_dm) num_test = num_dm;
  else srand(f->num_p);

  for (i=0; i<num_test; i++) {
    if (num_test == num_dm) j = i;
    else { j = rand(); j<<=31; j+=rand(); j%=(num_dm); }
    particle_r[i] = fast3tree_find_next_closest_distance(phasetree, res, 
							 f->particles[j].pos); //*pow(PARTICLE_MASS/f->particles[j].mass, 1.0/6.0);
  }
  target_r = find_median_r(particle_r, num_test, thresh);
  if (num_dm < f->num_p) fast3tree_rebuild(phasetree, f->num_p, f->particles);
  _find_subfofs_at_r(f, target_r);
}

void reassign_halo_particles(int64_t p_start, int64_t p_end) {
  int64_t last_halo, j;
  partition_sort_particles(p_start, p_end, copies, particle_halos);
  last_halo = particle_halos[p_start];
  halos[last_halo].p_start = p_start;
  for (j=p_start + 1; j<p_end; j++) {
    if (particle_halos[j] == last_halo) continue;
    halos[last_halo].num_p = j - halos[last_halo].p_start;
    last_halo = particle_halos[j];
    halos[last_halo].p_start = j;
  }
  halos[last_halo].num_p = j - halos[last_halo].p_start;
}

int could_be_poisson_or_force_res(struct halo *h1, struct halo *h2, int64_t *is_force_res) {
  float dx, r=0, v=0, mpe, mve; //, vt1, vt2;
  int64_t k;
  *is_force_res = 0;
  if (!h1->min_pos_err || !h1->min_vel_err) return 1;
  for (k=0; k<3; k++) {
    dx = h1->pos[k]-h2->pos[k];      r+=dx*dx;
    dx = h1->pos[k+3]-h2->pos[k+3];  v+=dx*dx;
  }
  mpe = h1->min_pos_err;
  mve = h1->min_vel_err;
  dx = (r / mpe + v / mve) / 2.0;
  if (!(dx > 100)) return 1;
  //if (h1->type != RTYPE_DM && h2->type != RTYPE_DM && !(dx>1000)) return 1;
  /*if (h1->type != RTYPE_DM && h2->type != RTYPE_DM) {
    float dist = _calc_halo_dist(h2, h1);
    if (h1->peak_density < 2.0*h2->av_density*exp(-dist*dist/2.0)) return 1;
    //printf("Marked halo %ld as poisson fluctuation (np=%"PRId64")\n",
    //h1-halos, h1->num_p);
    }*/

  r = sqrt(r);
  v = sqrt(v);
  if ((h1->r+h2->r > r) && (1.5*(h1->vrms+h2->vrms)) > v &&
      (r < 5*FORCE_RES)) {
    *is_force_res = 1;
    return 1;
  }
  return 0;
}

int64_t _find_biggest_parent(int64_t h_start, int64_t use_temporal_info,
			     int64_t growing) {
  int64_t i, j, max_i = h_start, num_m1, num_m2;
  float m1 = -1, max_vmax = 0, dx,ds, min_ds=0;
  if (growing) {
    assert(num_growing_halos);
    max_i = growing_halos[0]-halos;
  }
  for (i=h_start; i<num_halos; i++) {
    if (growing && !(halos[i].flags & GROWING_FLAG)) continue;
    if (halos[max_i].type == RTYPE_DM && halos[i].type != RTYPE_DM) continue;
    if (halos[i].vmax > halos[max_i].vmax) max_i = i;
    if (halos[i].vmax == halos[max_i].vmax &&
	halos[i].num_p > halos[max_i].num_p) max_i = i;
  }

  max_vmax = halos[max_i].vmax;
  if (use_temporal_info && TEMPORAL_HALO_FINDING && PARALLEL_IO) {
    for (i=h_start; i<num_halos; i++) {
      if (i==max_i) continue;
      if (halos[max_i].type == RTYPE_DM && halos[i].type != RTYPE_DM) continue;
      if (growing && !(halos[i].flags & GROWING_FLAG)) continue;
      if (max_vmax*0.6 < halos[i].vmax) {
	for (dx=0,j=0; j<3; j++) { ds=halos[i].pos[j]-halos[max_i].pos[j]; dx+=ds*ds; }
	if (!min_ds || min_ds > dx) min_ds = dx;
      }
    }
    min_ds = sqrt(min_ds)/3.0;

    for (i=h_start; i<num_halos; i++) {
      if (i==max_i) continue;
      if (growing && !(halos[i].flags & GROWING_FLAG)) continue;
      if (max_vmax*0.6 < halos[i].vmax) {
	if (m1 < 0)
	  m1 = find_previous_mass(halos+max_i, copies+halos[max_i].p_start,
				  &num_m1, min_ds);
	float m2 = find_previous_mass(halos+i, copies + halos[i].p_start,
				      &num_m2, min_ds);
	if (m1 && m2 && ((m2 > m1) || ((m2 == m1) && (num_m2 > num_m1)))) {
	  max_i = i;
	  m1 = m2;
	  num_m1 = num_m2;
	}
      }
    }
  }  
  return max_i;
}

void _fix_parents(int64_t h_start) {
  int64_t i, j, sub_of, num_m1, num_m2;
  float m1, m2, dx, ds;

  for (i=h_start; i<num_halos; i++) {
    extra_info[i].next_cochild = extra_info[i].prev_cochild = 
      extra_info[i].child = -1;
    if (extra_info[i].sub_of == i) extra_info[i].sub_of = -1;
  }

  if (TEMPORAL_HALO_FINDING && PARALLEL_IO) {
    for (i=h_start; i<num_halos; i++) {
      sub_of = extra_info[i].sub_of;
      if (sub_of == i) sub_of = extra_info[i].sub_of = -1;
      if (sub_of > -1 && halos[i].vmax > 0.6*halos[sub_of].vmax &&
	  (halos[i].type == RTYPE_DM || halos[sub_of].type != RTYPE_DM)) {
	for (dx=0,j=0; j<3; j++) { ds=halos[i].pos[j]-halos[sub_of].pos[j]; dx+=ds*ds; }
	dx = sqrt(dx)/3.0;
	m2 = find_previous_mass(halos+i, copies+halos[i].p_start, &num_m2, dx);
	m1 = find_previous_mass(halos+sub_of, copies+halos[sub_of].p_start,
				&num_m1, dx);
	if (m1 && m2 && ((m2 > m1) || ((m2 == m1) && (num_m2 > num_m1)))) {
	  extra_info[i].sub_of = extra_info[sub_of].sub_of;
	  extra_info[sub_of].sub_of = i;
	  if (extra_info[sub_of].sub_of == sub_of) 
	    extra_info[sub_of].sub_of = -1;
	  i--;
	}
      }
    }    
  }

  for (i=h_start; i<num_halos; i++) {
    extra_info[i].max_metric = 1e10;
    sub_of = extra_info[i].sub_of;
    if (sub_of > -1)
      extra_info[i].max_metric = _calc_halo_dist(halos+i, halos + sub_of);
    else continue;
    int64_t next_child = extra_info[sub_of].child;
    extra_info[i].next_cochild = next_child;
    extra_info[i].prev_cochild = -1;
    extra_info[sub_of].child = i;
    if (next_child > -1)
      extra_info[next_child].prev_cochild = i;
  }
}

void output_level(int64_t p_start, int64_t p_end, int64_t h_start, int64_t level)
{
  int64_t i;
  char buffer[1024];
  snprintf(buffer, 1024, "%s/levels_%f", OUTBASE, SCALE_NOW);
  FILE *output = check_fopen(buffer, "a");
  for (i=p_start; i<p_end; i++) {
    fprintf(output, "%f %f %f %f %f %f %"PRId64" %"PRId64" %"PRId64" %"PRId64" %"PRId32"\n",
	    copies[i].pos[0], copies[i].pos[1], copies[i].pos[2], 
	    copies[i].pos[3], copies[i].pos[4], copies[i].pos[5], 
	    p[copies[i].id].id, particle_halos[i], level, 
	    extra_info[particle_halos[i]].sub_of, copies[i].type);
  }
  fclose(output);

  snprintf(buffer, 1024, "%s/halos_%f.levels", OUTBASE, SCALE_NOW);
  output = check_fopen(buffer, "a");
  for (i=h_start; i<num_halos; i++) {
    fprintf(output, "%f %f %f %f %f %f %"PRId64" %"PRId64" %f %f %f %f %"PRId64" %"PRId64" %"PRId64"\n",
	    halos[i].pos[0], halos[i].pos[1], halos[i].pos[2], 
	    halos[i].bulkvel[0], halos[i].bulkvel[1], halos[i].bulkvel[2], 
	    halos[i].num_p, halos[i].num_child_particles, 
	    halos[i].r, halos[i].vrms, sqrt(halos[i].min_pos_err), 
	    sqrt(halos[i].min_vel_err), i, extra_info[i].sub_of, level);
  }
  fclose(output);
}

int sort_by_sub_of(const void *a, const void *b) {
  const struct halo *c = *((const struct halo **)a);
  const struct halo *d = *((const struct halo **)b);
  int64_t sc = extra_info[c-halos].sub_of;
  int64_t sd = extra_info[d-halos].sub_of;
  if (sc < sd) return -1;
  if (sd < sc) return 1;
  if (c->vmax < d->vmax) return -1;
  if (d->vmax < c->vmax) return 1;
  return 0;
}

void merge_galaxies(int64_t num_dm_halos) {
  qsort(growing_halos+num_dm_halos, num_growing_halos-num_dm_halos,
	sizeof(struct halo *), sort_by_sub_of);
  int64_t main_h = growing_halos[num_dm_halos]-halos;
  for (int64_t i=num_dm_halos+1; i<num_growing_halos; i++) {
    int64_t cur_h = growing_halos[i]-halos;
    assert(extra_info[main_h].sub_of <= extra_info[cur_h].sub_of);
    if (extra_info[cur_h].sub_of != extra_info[main_h].sub_of) {
      main_h = growing_halos[i]-halos;
      continue;
    }
    assert(main_h != cur_h);
    for (int64_t j=0; j<halos[cur_h].num_p; j++)
      particle_halos[halos[cur_h].p_start+j] = main_h;
    halos[cur_h].num_p = 0;
  }  
}

void _find_subs(struct fof *f, int64_t level) {
  int64_t f_start, f_end, h_start, i, j, f_index;
  int64_t p_start, max_i = 0, is_force_res;

  //Verify that the number of DM particles is >MIN_HALO_PARTICLES
  if (f->num_p < MIN_HALO_PARTICLES) return;
  /*if (!level) {
    for (i=0,j=0; i<f->num_p; i++) if (f->particles[i].type==RTYPE_DM) j++;
    if (j < MIN_HALO_PARTICLES) return;
    }*/

  //Find subFOFs
  p_start = f->particles - copies;
  f_index = f - subfofs;
  _find_subfofs_better2(f, FOF_FRACTION);
  f_start = num_subfofs;
  copy_fullfofs(&subfofs, &num_subfofs, &num_alloced_subfofs);
  f_end = num_subfofs;

  h_start = num_halos;
  for (i=f_start; i<f_end; i++)
    if (subfofs[i].num_p > MIN_HALO_PARTICLES)
      _find_subs(subfofs + i, level+1);

  //Convert particle positions back to normal
  if (level>0) f = subfofs + f_index;
  for (j=0; j<f->num_p; j++) particle_halos[p_start + j] = -1;
  for (i=h_start; i<num_halos; i++)
    for (j=0; j<halos[i].num_p; j++)
      particle_halos[halos[i].p_start + j] = i;

  for (j=0; j<f->num_p; j++) {
    struct particle *c = f->particles + j;
    memcpy(c->pos, p[c->id].pos, sizeof(float)*6);
  }

  if (h_start == num_halos) {
    add_new_halo(); //New seed halo
    halos[h_start].type = (level>0) ? f->particles[0].type : RTYPE_DM;
  }
  max_i = _find_biggest_parent(h_start, 1, 0);

  num_growing_halos=1;
  if (num_growing_halos >= num_alloc_gh) add_more_growing_halos();
  halos[max_i].flags |= GROWING_FLAG;
  extra_info[max_i].sub_of = -1;
  growing_halos[0] = halos + max_i;
  for (i=h_start; i<num_halos; i++) {
    if ((i==max_i) || !(halos[i].flags & GROWING_FLAG)) continue;
    if (could_be_poisson_or_force_res(halos+i, halos+max_i, &is_force_res)) {
      halos[i].flags -= (halos[i].flags & GROWING_FLAG);
      halos[i].flags |= GALAXY_INELIGIBLE_FLAG;
      extra_info[i].sub_of = max_i;
      if (!is_force_res) {
	for (j=0; j<halos[i].num_p; j++)
	  particle_halos[halos[i].p_start+j] = max_i;
	halos[i].num_p = 0;
      }
      continue;
    }
    if (num_growing_halos >= num_alloc_gh) add_more_growing_halos();
    growing_halos[num_growing_halos] = halos+i;
    num_growing_halos++;
  }

  if (num_growing_halos==1) {
    for (j=0; j<f->num_p; j++)
      if (particle_halos[p_start + j] < 0)
	particle_halos[p_start + j] = max_i;
  } else {
    build_subtree(growing_halos, num_growing_halos);
    for (j=0; j<f->num_p; j++) {
      if (particle_halos[p_start + j] < 0) {
	struct halo *h = find_best_halo(copies+p_start+j, halos+max_i);
	particle_halos[p_start + j] = h - halos;
	while (extra_info[h-halos].sub_of > -1) {
	  float max_metric = extra_info[h-halos].max_metric;
	  if (calc_particle_dist(h, copies+p_start+j) > max_metric) {
	    particle_halos[p_start + j] = extra_info[h-halos].sub_of;
	    h = halos + extra_info[h-halos].sub_of;
	  }
	  else break;
	}
      }
    }
  }

  if (TEMPORAL_HALO_FINDING && PARALLEL_IO) {
    for (i=h_start; i<num_halos; i++) {
      if (extra_info[i].ph < 0 || extra_info[i].sub_of < 0) continue;
      int64_t sub_of = extra_info[i].sub_of;
      if (extra_info[i].ph == extra_info[sub_of].ph &&
	  extra_info[i].max_metric < 1.5) {
	extra_info[i].ph = -1;
	for (j=halos[i].p_start; j<halos[i].num_p+halos[i].p_start; j++)
	  particle_halos[j] = sub_of;
	halos[i].num_p = halos[i].r = halos[i].vrms = 0;
      }
    }
  }
  reassign_halo_particles(p_start, p_start + f->num_p);
  calc_num_child_particles(h_start);
  for (i=0; i<num_growing_halos; i++) calc_basic_halo_props(growing_halos[i]);
  max_i = _find_biggest_parent(h_start, 0, 1);
  build_subtree(growing_halos, num_growing_halos);
  for (i=0; i<num_growing_halos; i++)
    extra_info[growing_halos[i]-halos].sub_of = 
      find_best_parent(growing_halos[i], halos+max_i) - halos;
  _fix_parents(h_start);
  calc_num_child_particles(h_start);
  for (i=0; i<num_growing_halos; i++) calc_basic_halo_props(growing_halos[i]);
  if (OUTPUT_LEVELS) output_level(p_start, p_start+f->num_p, h_start, level);


  num_subfofs = f_start;
  if (!level) {
    num_growing_halos = num_halos - h_start;
    while (num_alloc_gh < num_growing_halos) add_more_growing_halos();
    int64_t num_dm_halos = 0;
    for (i=h_start; i<num_halos; i++) {
      if (halos[i].type == RTYPE_DM) {
	growing_halos[num_dm_halos++] = halos+i;
	continue;
      }
      growing_halos[num_growing_halos-1-(i-num_dm_halos-h_start)] = halos+i;
    }
    
    for (i=0; i<num_growing_halos; i++)
      calc_basic_halo_props(growing_halos[i]);
    if (!num_dm_halos) {
      max_i = _find_biggest_parent(h_start, 0, 0);
      halos[max_i].type = RTYPE_DM;
      num_dm_halos = 1;
      assert(growing_halos[num_growing_halos-1-(max_i-h_start)] == halos+max_i);
      growing_halos[num_growing_halos-1-(max_i-h_start)] = growing_halos[0];
      growing_halos[0] = halos+max_i;
    }
    build_subtree(growing_halos, num_dm_halos);
    max_i = _find_biggest_parent(h_start, 0, 0);
    for (i=0; i<num_growing_halos; i++) {
      extra_info[growing_halos[i]-halos].sub_of = 
	find_best_parent(growing_halos[i], halos+max_i) - halos;
      /*      if (growing_halos[i] == halos+1210) {
	printf("(%f %f %f; %f %f %f; %f; %f)\n",
	       growing_halos[i]->pos[0], growing_halos[i]->pos[1], growing_halos[i]->pos[2],
	       growing_halos[i]->bulkvel[0],	       growing_halos[i]->bulkvel[1],	       growing_halos[i]->bulkvel[2],
	       growing_halos[i]->vmax, 	       growing_halos[i]->vmax_r);
	printf("(%f %f %f; %f %f %f; %f; %f)\n",
	       growing_halos[i]->pos[0], growing_halos[i]->pos[1], growing_halos[i]->pos[2],
	       growing_halos[i]->corevel[0],	       growing_halos[i]->corevel[1],	       growing_halos[i]->corevel[2],
	       growing_halos[i]->vmax, 	       growing_halos[i]->vmax_r);
      }
      if (dist(growing_halos[i]->pos, halos[1210].pos) < 0.03 &&
	  growing_halos[i]->type != RTYPE_DM) {
	printf("%ld; %f (%f %f %f; %f %f %f; %f; %f); %"PRId64"; %f\n",
	       growing_halos[i]-halos, calc_halo_dist(halos+1210, growing_halos[i]),
	       growing_halos[i]->pos[0], growing_halos[i]->pos[1], growing_halos[i]->pos[2],
	       growing_halos[i]->bulkvel[0],	       growing_halos[i]->bulkvel[1],	       growing_halos[i]->bulkvel[2],
	       growing_halos[i]->vmax, 	       growing_halos[i]->vmax_r,
	       extra_info[growing_halos[i]-halos].sub_of, 
	       calc_halo_dist(halos+extra_info[growing_halos[i]-halos].sub_of,
			      growing_halos[i]));
			      }*/
    }
    _fix_parents(h_start);

    merge_galaxies(num_dm_halos);
    reassign_halo_particles(p_start, p_start + f->num_p);
    for (i=0; i<num_growing_halos; i++) {
      calc_basic_halo_props(growing_halos[i]);
      convert_and_sort_core_particles(growing_halos[i], 
	    copies + growing_halos[i]->p_start, 0, NULL);
    }
    calc_num_child_particles(h_start);
    if (OUTPUT_LEVELS) output_level(p_start, p_start+f->num_p, h_start, -1);
  }
}

void find_subs(struct fof *f) {
  struct fof cf;
  int64_t i, h_start = num_halos;

  if (!phasetree) phasetree = fast3tree_init(0, NULL);
  if (!res) res = fast3tree_results_init();

  if (f->num_p > num_alloc_pc) alloc_particle_copies(f->num_p);
  if (!f->num_p) return;
  memcpy(copies, f->particles, sizeof(struct particle)*f->num_p);
  for (i=0; i<f->num_p; i++) copies[i].id = (f->particles-p)+i; //Hijack particle IDs
  cf = *f;
  cf.particles = copies;
  num_copies = f->num_p;

  if (LIGHTCONE) lightcone_set_scale(f->particles->pos);

  num_subfofs = 0;
  _find_subs(&cf, 0);
  num_subfofs = 0;
  for (i=0; i<f->num_p; i++) copies[i] = p[copies[i].id];
  calc_num_child_particles(h_start);
  for (i=h_start; i<num_halos; i++) calc_basic_halo_props(halos + i);
  for (i=h_start; i<num_halos; i++) calc_additional_halo_props(halos + i);

  memcpy(f->particles, copies, sizeof(struct particle)*f->num_p);
  for (i=h_start; i<num_halos; i++)
    halos[i].p_start += (f->particles - p);
}


void alloc_particle_copies(int64_t total_copies) {
  int64_t max_particle_r = MAX_PARTICLES_TO_SAMPLE;
  if (EXACT_LL_CALC) max_particle_r = total_copies;
  if (total_copies - num_alloc_pc < 1000) total_copies = num_alloc_pc + 1000;
  check_realloc_s(copies, sizeof(struct particle), total_copies);
  check_realloc_s(particle_halos, sizeof(int64_t), total_copies);
  if (max_particle_r > total_copies) max_particle_r = total_copies;
  check_realloc_s(particle_r, sizeof(float), max_particle_r);
  check_realloc_s(po, sizeof(struct potential), total_copies);
  num_alloc_pc = total_copies;
}

void free_particle_copies(void) {
  num_alloc_pc = 0;
  copies = check_realloc(copies, 0, "Freeing copies.");
  particle_halos = check_realloc(particle_halos, 0, "Freeing particle links.");
  particle_r = check_realloc(particle_r, 0, "Freeing particle radii.");
  po = check_realloc(po, 0, "Freeing potentials.");
  free_subtree();
}

void free_halos(void) {
  check_realloc_s(halos, 0, 0);
  num_halos = 0;
  if (num_alloc_gh) {
    check_realloc_s(growing_halos, 0, 0);
    num_alloc_gh = 0;
  }
  if (num_alloced_halo_ids) {
    check_realloc_s(halo_ids, 0, 0);
    num_alloced_halo_ids = 0;
  }
  check_realloc_s(extra_info, 0, 0);
}


int64_t rad_partition(float *rad, int64_t left, int64_t right, int64_t pivot_ind) {
  float pivot = rad[pivot_ind], tmp;
  int64_t si, i;
#define SWAP(a,b) { tmp = rad[a]; rad[a] = rad[b]; rad[b] = tmp; }
  SWAP(pivot_ind, right-1);
  si = right-2;
  for (i = left; i<si; i++) {
    if (rad[i] > pivot) { SWAP(i, si); si--; i--; }
  }
  if (rad[si] < pivot) si++;
  SWAP(right-1, si);
  return si;
#undef SWAP
}

inline float random_unit(void) {
  return(((float)(rand()%(RAND_MAX))/(float)(RAND_MAX)));
}

float find_median_r(float *rad, int64_t num_p, float frac) {
  int64_t pivot_index, k = num_p * frac;
  int64_t left = 0, right = num_p;
  assert(num_p>0);
  if (num_p < 2) return rad[0];
  while (1) {
    pivot_index = rad_partition(rad, left, right, 
				left + random_unit()*(right-left));
    if (k == pivot_index || (rad[left]==rad[right-1])) return rad[k];
    if (k < pivot_index) right = pivot_index;
    else left = pivot_index+1;
  }
}


void norm_sd(struct fof *f, float thresh) {
  double pos[6] = {0};
  double corr[6][6];
  double sig_x, sig_v, max_x=0, min_x=0, offset_x = 0;
  int64_t i,j,k,dm=1;
  double total_mass = 0;

  if (!f->num_p) return;

  for (i=0; i<f->num_p; i++)
    if (f->particles[i].type == RTYPE_DM) break;
  if (i==f->num_p) dm = 0;

  for (i=0; i<f->num_p; i++) {
    if (dm && f->particles[i].type != RTYPE_DM) continue;
    total_mass += f->particles[i].mass;
    for (j=0; j<6; j++) pos[j]+=f->particles[i].pos[j]*f->particles[i].mass;
  }
  if (!total_mass) return;

  for (j=0; j<6; j++) pos[j]/=total_mass;
  for (j=0; j<6; j++) for (k=0; k<6; k++) corr[j][k] = 0;

  for (i=0; i<f->num_p; i++) {
    for (j=0; j<6; j++) {
      f->particles[i].pos[j] -= pos[j];
      if (dm && f->particles[i].type != RTYPE_DM) continue;
      for (k=j; k<6; k++)
	corr[j][k]+=f->particles[i].pos[j]*f->particles[i].pos[k]*
	  f->particles[i].mass;
    }
  }

  for (j=0; j<6; j++)
    for (k=j; k<6; k++)
      corr[j][k]/=total_mass;

  calc_deviations(corr, &sig_x, &sig_v);

  if (f->num_p == num_copies) 
    sig_x *= INITIAL_METRIC_SCALING;
  if (!dm) sig_x /= NON_DM_METRIC_SCALING;

  if (!sig_x || !sig_v) return;

  for (i=0; i<f->num_p; i++) {
    for (j=0; j<6; j++)
      f->particles[i].pos[j] /= ((j < 3) ? sig_x : sig_v);
    if (f->particles[i].pos[0] > max_x) max_x = f->particles[i].pos[0];
    if (f->particles[i].pos[0] < min_x) min_x = f->particles[i].pos[0];
  }
  
  offset_x = 2.0*(max_x - min_x) + 1.0;
  /*for (i=0; i<f->num_p; i++)
    if (f->particles[i].type != RTYPE_DM)
    f->particles[i].pos[0] += offset_x;*/
  for (i=0; i<f->num_p; i++)
    f->particles[i].pos[0] += offset_x*f->particles[i].type;
}
