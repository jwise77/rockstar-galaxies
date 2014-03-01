#ifndef _PARENTS_C_
#define _PARENTS_C_

#ifndef GROUP_LIST
#error Must define GROUP_LIST!
#endif

#ifndef _FAST3TREE_C_
#error Must first include fast3tree.c!
#endif

#ifndef RADIUS_CONVERSION
#define RADIUS_CONVERSION 1.0e-3
#endif

#undef sort_halo_order
#define sort_halo_order _F3TN(FAST3TREE_PREFIX,sort_halo_order)
int sort_halo_order(const void *a, const void *b) {
  float c = GROUP_LIST[*((int64_t *)a)].vmax;
  float d = GROUP_LIST[*((int64_t *)b)].vmax;
  if (c > d) return -1;
  if (d > c) return 1;
  return 0;
}

#undef find_parents
#define find_parents _F3TN(FAST3TREE_PREFIX,find_parents)
void find_parents(int64_t ngroups) {
  int64_t i, j, *halo_order = NULL;
  struct fast3tree_results *nearest;
  struct fast3tree *halo_tree;
  FAST3TREE_TYPE *h1, *h2;
#ifndef NO_PERIODIC_BOUNDARIES
  float max_dist = BOX_SIZE/2.01;
#endif /*NO_PERIODIC_BOUNDARIES*/
  float range;

  halo_order = check_realloc(halo_order, sizeof(int64_t)*ngroups,
			     "Allocating halo order.");
  halo_tree = fast3tree_init(ngroups, GROUP_LIST);
#ifndef NO_PERIODIC_BOUNDARIES
  _fast3tree_set_minmax(halo_tree, 0, BOX_SIZE);
#endif /*NO_PERIODIC_BOUNDARIES*/
  for (i=0; i<ngroups; i++) {
    GROUP_LIST[i].parent = -1;
    halo_order[i] = i;
#ifdef TAG_MAJOR_MERGERS
    GROUP_LIST[i].mm = 0;
#endif /* TAG_MAJOR_MERGERS */
  }
  qsort(halo_order, ngroups, sizeof(int64_t), sort_halo_order);

  nearest = fast3tree_results_init();
  for (i=0; i<ngroups; i++) {
    h1 = &(GROUP_LIST[halo_order[i]]);
#ifdef IGNORE_INVALID_IDS
    if (h1->id < 0) continue;
#endif /* IGNORE_INVALID_IDS */

    range = h1->RADIUS*RADIUS_CONVERSION;
#ifndef NO_PERIODIC_BOUNDARIES
    if (max_dist < range) range = max_dist;
    fast3tree_find_sphere_periodic(halo_tree, nearest, h1->pos, range);
#else
    fast3tree_find_sphere(halo_tree, nearest, h1->pos, range);
#endif /* NO_PERIODIC_BOUNDARIES */

    for (j=0; j<nearest->num_points; j++) {
      h2 = nearest->points[j];
      if (h2->RADIUS < h1->RADIUS) {
	h2->parent = h1->id;
#ifdef TAG_MAJOR_MERGERS
	if (h2->vmax > 0.6*h1->vmax) h2->mm = h1->mm = 1;
#endif /*TAG_MAJOR_MERGERS */
      }
    }
  }
  fast3tree_results_free(nearest);
  fast3tree_free(&halo_tree);
  free(halo_order);
}

#endif /* _PARENTS_C_ */
