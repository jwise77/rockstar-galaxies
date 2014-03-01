#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <string.h>
#include <assert.h>
#include "../config_vars.h"
#include "../config.h"
#include "../rockstar.h"
#include "../io/meta_io.h"
#include "../io/io_bgc2.h"
#include "../check_syscalls.h"
#include "../io/bgc2.h"
#include "load_bgc2.h"
#include "../groupies.c"
#include "../halo.h"
#include "../potential.h"

float calc_distance2(float *p1, float *p2, float box_size) {
  double ds=0, dx;
  for (int64_t i=0; i<3; i++) {
    dx = fabs(p1[i]-p2[i]);
    if (dx > box_size/2.0) dx = box_size - dx;
    ds += dx*dx;
  }
  return ds;
}

int main(int argc, char **argv) {
  int64_t i,j,k, max_p=0;
  struct bgc2_header hdr;
  float bulkvel[3] = {0};
  double total_mass = 0;
  int64_t USE_COM = 0;

  if (argc < 9) {
    printf("Usage: %s mass_definition radius_fraction weighted_shapes force_resolution shape_iterations bound_particles use_com file1.bgc2 ...\n", argv[0]);
    printf("E.g.: %s vir 0.5 1 0.001 10 1 1 file1.bgc2 #calculates weighted shape for bound particles within 0.5 Rvir, assuming a force resolution of 0.001 Mpc/h, and does up to 10 iterations, using the halo center-of-mass within the chosen radius.\n", argv[0]);
    exit(1);
  }

  float rfrac = atof(argv[2]);
  WEIGHTED_SHAPES = atol(argv[3]);
  FORCE_RES = atof(argv[4]);
  SHAPE_ITERATIONS = atol(argv[5]);
  BOUND_PROPS = atol(argv[6]);
  USE_COM = atol(argv[7]);

  printf("#HaloID HaloParent Mass Radius Vmax X Y Z VX VY VZ Mass_within_%gxR%s A[x] A[y] A[z] b_to_a c_to_a\n", rfrac, argv[1]);
  printf("#Units: positions in Mpc/h (comoving)\n");
  printf("#Units: radii and other lengths in kpc/h (comoving)\n");
  printf("#Units: velocities in km/s (non-comoving)\n");
  printf("#HaloParent: the ID of a larger halo containing this halo, if any; -1 otherwise\n");
  printf("#Shapes calculated within %g x R%s\n", rfrac, argv[1]);
  printf("#Shapes are %s\n", (WEIGHTED_SHAPES ? "weighted" : "unweighted"));
  if (BOUND_PROPS) printf("#Only bound particles are used.\n");
  else printf("#Both bound and unbound particles are used.\n");
  printf("#Assumed force resolution (min. radius): %f Mpc/h\n", FORCE_RES);
  printf("#Maximum number of shape iterations: %"PRId64"\n", SHAPE_ITERATIONS);
  printf("#Using halo %s\n", (USE_COM ? "center of mass" : "density peak"));
  printf("#b_to_a, c_to_a: Ratio of second and third largest shape ellipsoid axes (B and C) to largest shape ellipsoid axis (A) (dimensionless).\n#  Shapes are determined by the method in Allgood et al. (2006).\n#A[x],A[y],A[z]: Largest shape ellipsoid axis (kpc/h comoving)\n#All particles are considered (bound and unbound) so these measurements\n#  may differ from the Rockstar catalogs, especially for subhalos.\n");
  for (i=8; i<argc; i++) {
    num_groups = num_parts = 0;
    load_bgc2(argv[i], &hdr, &grps, &num_groups, &parts, &num_parts);
    if (i==8) {
      printf("#Box size: %g Mpc/h\n", hdr.box_size);
      printf("#Particle Mass: %g Msun/h\n", hdr.part_mass);
      printf("#h: %g; Om: %g; Ol: %g\n", hdr.Hubble0, hdr.Omega0, hdr.OmegaLambda);
      printf("#Redshift: %g\n", hdr.redshift);
    }
    int64_t p_start = 0;
    SCALE_NOW = 1.0/(hdr.redshift+1.0);
    h0 = hdr.Hubble0;
    Om = hdr.Omega0;
    Ol = hdr.OmegaLambda;
    PARTICLE_MASS = hdr.part_mass;
    MASS_DEFINITION = argv[1];
    calc_mass_definition();
    double dens_thresh = particle_thresh_dens[0]*(4.0*M_PI/3.0)*hdr.part_mass;
    for (j=0; j<hdr.ngroups; j++) {
      struct halo h;
      memset(&h, 0, sizeof(struct halo));
      check_realloc_var(po, sizeof(struct potential), max_p, grps[j].npart);
      for (k=0; k<grps[j].npart; k++) {
	memcpy(po[k].pos, parts[p_start+k].pos, sizeof(float)*3);
	memcpy(po[k].pos+3, parts[p_start+k].vel, sizeof(float)*3);
	po[k].mass = parts[p_start+k].mass;
	for (int64_t l=0; l<3; l++) bulkvel[l] += po[k].pos[l+3]*po[k].mass;
	total_mass += po[k].mass;
	po[k].pe = po[k].ke = 0;
	po[k].flags = 0;
	po[k].r2 = calc_distance2(grps[j].pos, po[k].pos, hdr.box_size);
      }
      memcpy(h.pos, grps[j].pos, sizeof(float)*3);
      memcpy(h.pos+3, grps[j].vel, sizeof(float)*3);
      if (total_mass > 0)
 	for (int64_t l=0; l<3; l++) bulkvel[l] /= total_mass;
      if (BOUND_PROPS) {
 	compute_potential(po, grps[j].npart);
 	compute_kinetic_energy(po, grps[j].npart, bulkvel, h.pos);
      }
      qsort(po, grps[j].npart, sizeof(struct potential), dist_compare);

      h.r = 0;
      int64_t np = 0;
      for (k=0; k<grps[j].npart; k++) {
	total_mass += po[k].mass;
	double r = sqrt(po[k].r2);
	if (total_mass/(r*r*r) >= dens_thresh) { np = k; }
      }
      k = np;
      if (k>=0) h.r = po[k].r2;
      double pos[3] = {0};
      total_mass = 0;
      for (k=0; k<grps[j].npart; k++) {
	total_mass += po[k].mass;
	for (int64_t l=0; l<3; l++) pos[l] += po[k].pos[l]*po[k].mass;
	if (po[k].r2 > rfrac*rfrac*h.r) break;
      }

      if (USE_COM && total_mass>0)
 	for (int64_t l=0; l<3; l++) h.pos[l] = pos[l]/total_mass;
      h.r = rfrac*sqrt(h.r)*1e3;
      h.m = total_mass;
      calc_shape(&h, k, BOUND_PROPS);
      printf("%"PRId64" %"PRId64" %e %f %f %f %f %f %f %f %f %e %f %f %f %f %f\n",
	     grps[j].id, grps[j].parent_id, grps[j].mass, 1e3*grps[j].radius, grps[j].vmax, grps[j].pos[0], grps[j].pos[1], grps[j].pos[2], grps[j].vel[0], grps[j].vel[1], grps[j].vel[2], k*PARTICLE_MASS, h.A[0], h.A[1], h.A[2], h.b_to_a, h.c_to_a);
      p_start += grps[j].npart;
    }
  }
  return 0;
}
