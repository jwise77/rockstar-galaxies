#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include "config_vars.h"
#include "subhalo_metric.h"
#include "check_syscalls.h"
#include "groupies.h"

#define FAST3TREE_TYPE struct halo_metric
#define FAST3TREE_DIM 4
#define FAST3TREE_PREFIX SUBHALO_METRIC
#define POINTS_PER_LEAF 10
#include "fast3tree.c"

#define INV_RADIUS_WEIGHTING -0.2 //Forces divisions along the radius dimension

struct fast3tree *subtree = NULL;
struct halo_metric *sub_metric = NULL;
struct fast3tree_results *subtree_res = NULL;
int64_t alloced_metrics = 0;

void build_subtree(struct halo **subs, int64_t num_subs) {
  int64_t i;
  if (num_subs > alloced_metrics) {
    alloced_metrics = num_subs;
    sub_metric = check_realloc(sub_metric, sizeof(struct halo_metric)*num_subs,
			       "Allocating room for halo metric tree.");
  }
  for (i=0; i<num_subs; i++) {
    memcpy(sub_metric[i].pos, subs[i]->pos, sizeof(float)*3);
    sub_metric[i].pos[3] = subs[i]->r*(1.0/INV_RADIUS_WEIGHTING);
    sub_metric[i].target = subs[i];
  }
  if (subtree) fast3tree_rebuild(subtree, num_subs, sub_metric);
  else {
    subtree = fast3tree_init(num_subs, sub_metric);

    subtree_res = fast3tree_results_init();
  }
}

void free_subtree(void) {
  fast3tree_free(&subtree);
  alloced_metrics = 0;
  sub_metric = check_realloc(sub_metric, 0, "Freeing halo metric tree.\n");
}

float calc_particle_dist(struct halo *h, struct particle *part) {
  float dx, r2=0, v2=0; //,dv;
  int64_t i; //,j;
  float *vel = h->bulkvel;
  if (h->vrms <= 0 || h->r <= 0 || h->num_p <= 0) return 1e20;
  if (part->type == RTYPE_DM && h->type != RTYPE_DM) return 1e20;
  if (part->type != RTYPE_DM && h->type == RTYPE_DM &&
      h->flags & GALAXY_INELIGIBLE_FLAG) return 1e20;
  //if (h->type == RTYPE_DM) {
    //if (part->type != RTYPE_DM) vel = h->corevel;
    for (i=0; i<3; i++) { dx = h->pos[i]-part->pos[i]; r2+=dx*dx; }
    for (; i<6; i++)  { dx = vel[i-3]-part->pos[i]; v2+=dx*dx; }
    if (h->type == RTYPE_DM && part->type != RTYPE_DM) //Downweight velocity differences for gas + stars
      v2 /= (NON_DM_METRIC_SCALING*NON_DM_METRIC_SCALING);
    return sqrt((r2 / (h->r*h->r)) + v2 / (h->vrms*h->vrms));
    //}

    /*  struct extra_halo_info *ei = extra_info + (h-halos);
  for (i=0; i<3; i++) {
    for (dx=0,dv=0,j=0; j<3; j++) {
      dx += ei->x_orth_matrix[i][j]*(part->pos[j]-h->pos[j]);
      dv += ei->v_orth_matrix[i][j]*(part->pos[j+3]-h->pos[j+3]);
    }
    r2 += dx*dx/ei->x_eig[i];
    v2 += dv*dv/ei->v_eig[i];
  }
  v2 /= (NON_DM_METRIC_SCALING*NON_DM_METRIC_SCALING);
  return sqrt(r2 + v2);
    */
}

float _calc_halo_dist(struct halo *h1, struct halo *h2) {
  struct particle tmp;
  memcpy(tmp.pos, h2->pos, sizeof(float)*3);
  memcpy(tmp.pos+3, h2->bulkvel, sizeof(float)*3);
  tmp.type = h2->type;
  return calc_particle_dist(h1, &tmp);
}

float calc_halo_dist(struct halo *h1, struct halo *h2) {
  if (h2->r > h1->r*0.99999 && 
      !(h1->type == RTYPE_DM && h2->type != RTYPE_DM)) return 1e20;
  if ((h1->type == RTYPE_DM && h2->type != RTYPE_DM) && (h1->vmax < 0.2*h2->vmax)) return 1e20;
  //if (h2->r > h1->r*2.0) return 1e20;
  return (_calc_halo_dist(h1, h2));
}

int64_t node_could_be_better(struct tree3_node *n, struct particle *part,
			     float best_metric) {
  float max_r = n->min[3]*best_metric*INV_RADIUS_WEIGHTING;
  int64_t i;
  for (i=0; i<3; i++)
    if ((part->pos[i]+max_r < n->min[i]) || (part->pos[i]-max_r > n->max[i]))
      return 0;
  return 1;
}

struct halo *_find_best_halo(struct tree3_node *n, struct particle *part,
			     float *best_metric, struct halo *best_halo) {
  int64_t i;
  float metric;
  if (n->div_dim < 0) { //At leaf node
    for (i=0; i<n->num_points; i++) {
      metric = calc_particle_dist(n->points[i].target, part);
      if (metric < *best_metric) {
	*best_metric = metric;
	best_halo = n->points[i].target;
      }
    }
  } else {
    if (node_could_be_better(n->left, part, *best_metric))
      best_halo = _find_best_halo(n->left, part, best_metric, best_halo);
    if (node_could_be_better(n->right, part, *best_metric))
      best_halo = _find_best_halo(n->right, part, best_metric, best_halo);
  }
  return (best_halo);
}


double calc_expected_density(struct halo *h, struct particle *part) {
  double dx, r2 = 0, v2=0, r, prof;
  int64_t i;
  if (h->vrms <= 0 || h->r <= 0) return 1e9;
  for (i=0; i<6; i++) { dx = h->pos[i]-part->pos[i];
    if (i<3) r2+=dx*dx; else v2+=dx*dx; }
  r = sqrt(r2);
  prof = 1.0+10.0*r/h->r;
  return (exp(-0.5*(v2 / (h->vrms*h->vrms))) / (h->vrms*h->vrms*h->vrms) *
	  h->m/(h->r*h->r*r*(prof*prof)));
}

struct halo *alt_find_best_halo(struct particle *part, struct halo *best_halo) {
  double best_metric = calc_expected_density(best_halo, part);
  double density;
  int64_t i;
  for (i=0; i<subtree->root->num_points; i++) {
    density = calc_expected_density(subtree->root->points[i].target, part);
    if (density > best_metric) {
      best_metric = density;
      best_halo = subtree->root->points[i].target;
    }
  }
  return (best_halo);
}


struct halo *find_best_halo(struct particle *part, struct halo *best_halo) {
  float best_metric = calc_particle_dist(best_halo, part);
  if (ALT_NFW_METRIC) return (alt_find_best_halo(part, best_halo));
 return (_find_best_halo(subtree->root, part, &best_metric, best_halo));
}


int64_t node_could_be_better_parent(struct tree3_node *n, struct halo *h,
				    float best_metric) {
  float r = n->min[3]*INV_RADIUS_WEIGHTING;
  float max_r = r*best_metric;
  int64_t i;
  if ((r < h->r*1.00001) && (h->type==RTYPE_DM)) return 0;
  for (i=0; i<3; i++)
    if ((h->pos[i]+max_r < n->min[i]) || (h->pos[i]-max_r > n->max[i]))
      return 0;
  return 1;
}

struct halo *_find_best_parent(struct tree3_node *n, struct halo *h,
			       float *best_metric, struct halo *best_halo) {
  int64_t i;
  float metric;
  if (n->div_dim < 0) { //At leaf node
    for (i=0; i<n->num_points; i++) {
      if (h->type == RTYPE_DM && 
	  n->points[i].target->type != RTYPE_DM) continue;
      metric = calc_halo_dist(n->points[i].target, h);
      if (metric < *best_metric) {
	*best_metric = metric;
	best_halo = n->points[i].target;
      }
    }
  } else {
    if (node_could_be_better_parent(n->left, h, *best_metric))
      best_halo = _find_best_parent(n->left, h, best_metric, best_halo);
    if (node_could_be_better_parent(n->right, h, *best_metric))
      best_halo = _find_best_parent(n->right, h, best_metric, best_halo);
  }
  return (best_halo);
}


struct halo *find_best_parent(struct halo *h, struct halo *biggest_halo) {
  //Need to find closest larger halo
  if (h==biggest_halo) return h;
  float best_metric = calc_halo_dist(biggest_halo, h);
  return _find_best_parent(subtree->root, h, &best_metric, biggest_halo);
}


struct halo_metric **find_children(struct halo *h, struct halo *parent, 
				   float r, int64_t *num_children) {
  int64_t i, j, swap;
  float dx, ds, t_r;
  float pos[4];
  memcpy(pos, h->pos, sizeof(float)*3);
  pos[3] = h->r*(1.0/INV_RADIUS_WEIGHTING);
  t_r = r*sqrt(1.0 + 1.0/(INV_RADIUS_WEIGHTING*INV_RADIUS_WEIGHTING));
  fast3tree_find_sphere(subtree, subtree_res, pos, t_r);
  for (i=0; i<subtree_res->num_points; i++) {
    swap = 0;
    if (subtree_res->points[i]->target->r >= 
	h->r) swap = 1;
    for (ds=0,j=0; j<3; j++) {
      dx = subtree_res->points[i]->pos[j]-h->pos[j];
      ds += dx*dx;
    }
    if (ds >= r*r) swap = 1;
    if (parent && (calc_halo_dist(h, subtree_res->points[i]->target) >
		   calc_halo_dist(parent, subtree_res->points[i]->target)))
      swap = 1;
    if (swap) {
      subtree_res->points[i] = subtree_res->points[subtree_res->num_points-1];
      subtree_res->num_points--;
      i--;
    }
  }
  *num_children = subtree_res->num_points;
  return subtree_res->points;
}
