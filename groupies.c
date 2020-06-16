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
#define FAST3TREE_EXTRA_INFO float ll, density;
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


void link_particle_types(struct fof *f) {
  int64_t i;
  assert(f->num_p == num_copies);
  init_particle_smallfofs(f->num_p, f->particles);
  int64_t idx[NUM_RTYPES];
  for (i=0; i<NUM_RTYPES; i++) idx[i] = -1;
  for (i=0; i<f->num_p; i++) {
    int64_t type = f->particles[i].type;
    if (idx[type] < 0) { idx[type] = i; continue; }
    struct particle *tp[2] = {f->particles+idx[type], f->particles+i};
    link_particle_to_fof(f->particles+idx[type], 2, tp);
  }
  build_fullfofs();
}

void _find_subfofs_at_r(struct fof *f, float target_r) {
  int64_t i;
  init_particle_smallfofs(f->num_p, f->particles);
  
  for (i=0; i<f->num_p; i++) {
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
  norm_sd(f, thresh, NULL, NULL);
  //norm_sd_bary(f);
  int64_t num_dm = 0; //f->num_p;
  //int64_t num_dm = separate_dm(f);
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

int sort_by_ll(const void *a, const void *b) {
  const int64_t *c = a;
  const int64_t *d = b;
  if (phasetree->root[*c].ll < phasetree->root[*d].ll) return -1;
  if (phasetree->root[*c].ll > phasetree->root[*d].ll) return 1;
  return 0;
}


int64_t poisson_merge(float ll, float ll_saddle, int64_t np) {
  if (np < MIN_HALO_PARTICLES) return 1;
  //Require ~5sig:
  float rel_err = sqrt(1.0/(np-1.0));
  if (ll*(1.0+GALAXY_POISSON_SIGMA*rel_err) > ll_saddle) return 1;
  return 0;
}

int64_t join_subfofs_if_needed(int64_t sf1, int64_t sf2, struct tree3_node *n, int64_t *gp, float *lls,
			       struct tree3_node **rn, FILE *distances) {
  if (sf1 == sf2) return sf1;
  int64_t tmp = sf1;
  double dx=0, dv=0;
  if (gp[sf2] > gp[sf1]) { sf1 = sf2; sf2 = tmp; }
  if (distances) {
    for (int64_t k=0; k<6; k++) {
      double ds = 0.5*(rn[sf1]->max[k] + rn[sf1]->min[k] - rn[sf2]->max[k] - rn[sf2]->min[k]);
      if (k<3) dx += ds*ds;
      else dv += ds*ds;
    }
    fprintf(distances, "%"PRId64" %"PRId64" %"PRId64" %"PRId64" %e %e %e %e\n", sf1, gp[sf1], sf2, gp[sf2], sqrt(dx), sqrt(dv), n->ll, n->density);
  }
  
  if (poisson_merge(lls[sf2], n->ll, gp[sf2]) || n->density > DENSITY_MERGE_THRESH) {
    if (gp[sf2]>100 && distances) {
      fprintf(distances, "#Merged sf%"PRId64"(%"PRId64") and sf%"PRId64"(%"PRId64"); ll %e & %e; dx=%f; dv=%f; dens=%e\n", sf1, gp[sf1], sf2, gp[sf2], lls[sf2], n->ll, dx, dv, n->density);
    }
    merge_smallfofs(smallfofs+sf1, smallfofs+sf2);
    lls[sf2] = lls[sf1];
    gp[sf1] += gp[sf2];
    gp[sf2] = 0;
  }
  return sf1;
}

void _find_subfofs_better3(struct fof *f) {
  int64_t i, j, k;
  static int64_t sf_offset = 0;
  norm_sd_bary(f);
  //float v_scaling = sqrt(SCALE_NOW)*dynamical_time/NON_DM_METRIC_SCALING;
  fast3tree_rebuild(phasetree, f->num_p, f->particles);
  for (i=0; i<phasetree->num_nodes; i++) {
    double v = 1.0;
    struct tree3_node *n = phasetree->root + i;
    for (j=0; j<6; j++) v *= n->max[j]-n->min[j];
    n->ll = 1e20;
    n->density = 0;
    if (n->num_points && v>0) n->ll = pow(v/n->num_points, 1.0/6.0);
    else if (n->parent) n->ll = n->parent->ll;
    for (j=0; j<n->num_points; j++) n->density += n->points[j].mass;
    if (n->num_points && v>0) n->density /= v;
    else n->density = 0;
  }

  int64_t *process_order = NULL;
  check_realloc_s(process_order, sizeof(int64_t), phasetree->num_nodes);
  int64_t num_leaf_nodes = 0;
  for (i=0; i<phasetree->num_nodes; i++) {
    if (phasetree->root[i].div_dim >= 0) continue;
    process_order[num_leaf_nodes] = i;
    num_leaf_nodes++;
  }
  qsort(process_order, num_leaf_nodes, sizeof(int64_t), sort_by_ll);
  init_particle_smallfofs(f->num_p, f->particles);

  float *lls = NULL;
  check_realloc_s(lls, sizeof(float), num_leaf_nodes);
  int64_t *gp = NULL;
  check_realloc_s(gp, sizeof(int64_t), num_leaf_nodes);
  struct tree3_node **rn = NULL;
  check_realloc_s(rn, sizeof(struct tree3_node *), num_leaf_nodes);
  FILE *distances = NULL;
  if (OUTPUT_LEVELS) {
    distances = check_fopen("fof_distances.txt", "a");
    fprintf(distances, "#FOF1 NP1 FOF2 NP2 dx dv ll\n");
  }
  for (i=0; i<num_leaf_nodes; i++) {
    //fprintf(stderr, ".");
    struct tree3_node *n = phasetree->root +  process_order[i];
    if (n->ll > 1e10 || !n->num_points) continue;
    //Need to sort by roots, also put -1's at end
    partition_sort_particles(0, n->num_points, n->points, particle_smallfofs + (n->points - f->particles));
    int64_t j_start = 0;
    int64_t main_halo = -1;
    for (; j_start<n->num_points; j_start++) {
      int64_t sf = SMALLFOF_OF(n->points+j_start);
      if (sf > -1) {
	_collapse_smallfof(smallfofs+sf);
	sf = smallfofs[sf].root;
	if (main_halo < 0 || gp[sf] > gp[main_halo])
	  main_halo = sf;
      }
    }
    if (main_halo < 0) {
      main_halo = add_new_smallfof();
      lls[main_halo] = n->ll;
      gp[main_halo] = 0;
      rn[main_halo] = n;
    }

    for (j=j_start+1; j<n->num_points; j++) {
      int64_t sf = SMALLFOF_OF(n->points+j);
      if (sf < 0) {
	sf = SMALLFOF_OF(n->points+j) = main_halo;
	gp[sf]++;
      }
      else {
	_collapse_smallfof(smallfofs+sf);
	sf = smallfofs[sf].root;
	main_halo = join_subfofs_if_needed(main_halo, sf, n, gp, lls, rn, distances);
      }
    }

    float ll = n->ll;
    ll *= GALAXY_LINKING_LENGTH;
    for (j=0; j<n->num_points; j++) {
      //fast3tree_find_sphere_marked(phasetree, res, n->points[j].pos, ll, 0, 1);
      fast3tree_find_sphere(phasetree, res, n->points[j].pos, ll);
      int64_t sf = main_halo;
      for (k=0; k<res->num_points; k++) {
	int64_t sf2 = SMALLFOF_OF(res->points[k]);
	if (sf2 < 0) {
	  SMALLFOF_OF(res->points[k]) = sf;
	  gp[sf]++;
	} else {
	  _collapse_smallfof(smallfofs+sf2);
	  sf2 = smallfofs[sf2].root;
	  sf = main_halo = join_subfofs_if_needed(sf, sf2, n, gp, lls, rn, distances);
	}
      }
    }
  }

  if (OUTPUT_LEVELS) {
    FILE *out = check_fopen("stars_levels.txt", "a");
    fprintf(out, "#X Y Z VX VY VZ LL Density SmallFoF Mass\n");
    int64_t max_sf = -1;
    for (i=0; i<phasetree->num_nodes; i++) {
      struct tree3_node *n = phasetree->root + i;
      if (n->div_dim >= 0) continue;
      for (j=0; j<n->num_points; j++) {
	int64_t sf = SMALLFOF_OF(n->points + j);
	if (sf>0) sf = smallfofs[sf].root;
	if (sf > max_sf) max_sf = sf;
	fprintf(out, "%f %f %f %f %f %f %e %e %"PRId64" %e\n", n->points[j].pos[0], n->points[j].pos[1], n->points[j].pos[2], n->points[j].pos[3], n->points[j].pos[4], n->points[j].pos[5], n->ll, n->density, sf, n->points[j].mass);
      }
    }
    fclose(out);
    sf_offset += max_sf + 1;
  }

  if (distances) fclose(distances);
  free(lls);
  free(gp);
  free(rn);
  free(process_order);
  build_fullfofs();
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
    assert(sub_of != i);
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
  fprintf(output, "#X Y Z VX VY VZ ID Halo Level Sub_Of Type\n");
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
  fprintf(output, "#X Y Z VX VY VZ NP NCP R VRMS MPErr MVerr HID Sub_Of M Vmax Type Level\n");
  for (i=h_start; i<num_halos; i++) {
    fprintf(output, "%f %f %f %f %f %f %"PRId64" %"PRId64" %f %f %f %f %"PRId64" %"PRId64" %e %e %"PRId32" %"PRId64"\n",
	    halos[i].pos[0], halos[i].pos[1], halos[i].pos[2], 
	    halos[i].bulkvel[0], halos[i].bulkvel[1], halos[i].bulkvel[2], 
	    halos[i].num_p, halos[i].num_child_particles, 
	    halos[i].r, halos[i].vrms, sqrt(halos[i].min_pos_err), 
	    sqrt(halos[i].min_vel_err), i, extra_info[i].sub_of,
	    halos[i].m, halos[i].vmax, halos[i].type,
	    level);
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
  int64_t p_start, max_i = 0, is_force_res, do_higher_levels=1;

  //Verify that the number of DM particles is >MIN_HALO_PARTICLES
  if (f->num_p < MIN_HALO_PARTICLES) return;
  if (f->particles[0].type == RTYPE_STAR) {
    if (level==2) do_higher_levels = 0;
  }

  //Find subFOFs
  p_start = f->particles - copies;
  f_index = f - subfofs;
  if (do_higher_levels) {
    if (level) {
      if (f->particles[0].type != RTYPE_STAR)
	_find_subfofs_better2(f, FOF_FRACTION);
      else
	_find_subfofs_better3(f);
    }
    else link_particle_types(f);
  }
  f_start = num_subfofs;
  if (do_higher_levels) 
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
    
    if (1) {
      int64_t pass = 0;
      for (pass=0; pass<1; pass++) {
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
	/*	for (i=0; i<num_growing_halos; i++) {
	  extra_info[growing_halos[i]-halos].sub_of = 
	    find_best_parent(growing_halos[i], halos+max_i) - halos;
	    }*/
	for (i=h_start; i<num_halos; i++) {
	  extra_info[i].sub_of = 
	    find_best_parent(halos+i, halos+max_i) - halos;
	}
	
	if (OUTPUT_LEVELS) {
	  FILE *out = check_fopen("galaxy_metrics.dat", "a");
	  fprintf(out, "#ID X Y Z VX VY VZ M V HID HX HY HZ HVX HVY HVZ HM HV Dist Dx/r Dv/vrms Sub_Of\n");
	  for (i=h_start; i<num_halos; i++) {
	    if (halos[i].type != RTYPE_STAR) continue;
	    for (j=0; j<num_growing_halos; j++) {
	      double dist = calc_halo_dist(growing_halos[j], halos+i);
	      if (dist < 1e10) {
		int64_t k;
		double dx=0, dv=0;
		for (k=0; k<3; k++) {
		  double ds = halos[i].pos[k]-growing_halos[j]->pos[k];
		  dx+=ds*ds;
		  ds = halos[i].pos[k+3]-growing_halos[j]->pos[k+3];
		  dv += ds*ds;
		}
		dx = sqrt(dx);
		dv = sqrt(dv);
		if (growing_halos[j]->r>0) dx /= growing_halos[j]->r;
		if (growing_halos[j]->vrms>0) dv /= growing_halos[j]->vrms;
		fprintf(out, "%"PRId64" %f %f %f %f %f %f %e %f %"PRId64" %f %f %f %f %f %f %e %f %e %f %f %"PRId64"\n", halos[i].id, halos[i].pos[0], halos[i].pos[1], halos[i].pos[2], halos[i].pos[3], halos[i].pos[4], halos[i].pos[5], halos[i].m, halos[i].vmax, growing_halos[j]->id, growing_halos[j]->pos[0], growing_halos[j]->pos[1], growing_halos[j]->pos[2], growing_halos[j]->pos[3], growing_halos[j]->pos[4], growing_halos[j]->pos[5], growing_halos[j]->m, growing_halos[j]->vmax, dist, dx, dv, extra_info[i].sub_of);
	      }
	    }
	  }
	  fclose(out);
	}
	_fix_parents(h_start);
      }
    }

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


void norm_sd_bary(struct fof *f) {
  float v_scaling = sqrt(SCALE_NOW)*dynamical_time/NON_DM_METRIC_SCALING/10.0;
  int64_t i, j;
  for (i=0; i<f->num_p; i++) {
    for (j=3; j<6; j++)
      f->particles[i].pos[j] *= v_scaling;
  }
}


void norm_sd(struct fof *f, float thresh, double *axis_x, double *axis_v) {
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
    for (j=0; j<6; j++)
      f->particles[i].pos[j] -= pos[j];
    if (dm && f->particles[i].type != RTYPE_DM) continue;
    for (j=0; j<6; j++)
      for (k=j; k<6; k++)
	corr[j][k]+=f->particles[i].pos[j]*f->particles[i].pos[k]*
	  f->particles[i].mass;
  }

  for (j=0; j<6; j++)
    for (k=j; k<6; k++)
      corr[j][k]/=total_mass;

  calc_deviations(corr, &sig_x, &sig_v, axis_x, axis_v);

  if (f->num_p == num_copies) 
    sig_x *= INITIAL_METRIC_SCALING;
  //if (!dm) sig_x /= NON_DM_METRIC_SCALING;
  //  float sig_x_bary = sig_x / 3.0;

  float v_scaling = 1.0/(sqrt(SCALE_NOW)*dynamical_time);
  sig_x = 1;
  if (f->num_p == num_copies)
    sig_v = sig_x*v_scaling;
  else {
    sig_x = sig_v = 1;
  }

  float sig_v_bary = sig_v * NON_DM_METRIC_SCALING;

  if (!sig_x || !sig_v) return;

  for (i=0; i<f->num_p; i++) {
    //float t_sig_v = (f->particles[i].type == RTYPE_DM) ? sig_v : sig_v_bary;
    float t_sig_v = sig_v_bary;
    for (j=0; j<6; j++)
      f->particles[i].pos[j] /= ((j < 3) ? sig_x : t_sig_v);
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
