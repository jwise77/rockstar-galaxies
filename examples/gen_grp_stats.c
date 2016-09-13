#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <math.h>
#include "../halo.h"
#include "../particle.h"
#include "../potential.h"
#include "../check_syscalls.h"
#include "load_full_particles.h"
#include <assert.h>
#include <sys/stat.h>
#include "../config_vars.h"
#include "../io/io_nchilada.h"
#include "../io/io_util.h"

#ifndef CALC_POTENTIALS
#error Need to define CALC_POTENTIALS when compiling this file!
#endif /* !def CALC_POTENTIALS */

extern int64_t NCHILADA_SWAP;
extern int64_t NCHILADA_ID_BYTES;

struct halo *h = NULL;
int64_t *parents = NULL;
struct full_particle *p = NULL;
float bounds[6] = {0};
int64_t num_h=0, num_p=0;
struct potential *po = NULL;
int64_t max_num_po = 0;
int64_t *grp_ids = NULL;
float *grp_masses = NULL;
int64_t num_grp_ids = 0;
FILE *stats, *grp, *ids;


void check_grp_ids_space(int64_t num_requested) {
  int64_t i;
  if (num_requested < num_grp_ids) return;
  check_realloc_s(grp_ids, sizeof(int64_t), num_requested);
  check_realloc_s(grp_masses, sizeof(float), num_requested);
  for (i=num_grp_ids; i<num_requested; i++) { grp_ids[i] = -1; grp_masses[i] = -1; }
  num_grp_ids = num_requested;
}

void calc_potentials(double dm_mass) {
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
      for (int64_t k=0; k<3; k++) {
	double dx = the_h->pos[k]-po[j].pos[k];
	if (dx > BOX_SIZE / 2.0) po[j].pos[k] += BOX_SIZE;
	if (dx < -BOX_SIZE / 2.0) po[j].pos[k] -= BOX_SIZE;
      }
      po[j].pe = po[j].ke = 0;
      po[j].mass = p[the_h->p_start+j].mass;
      po[j].energy = p[the_h->p_start+j].energy;
      po[j].hid = p[the_h->p_start+j].ehid;
      po[j].type = p[the_h->p_start+j].type;
    }
    compute_kinetic_energy(po, the_h->num_p, the_h->pos+3, the_h->pos);
    compute_potential(po, the_h->num_p);
    count=0;

    int64_t n[3]={0};
    double m[3]={0};

    for (j=0; j<the_h->num_p; j++) {
      if (po[j].pe < po[j].ke) continue;
      double r=0;
      for (int64_t k=0; k<3; k++) {
	double dx = the_h->pos[k]-po[j].pos[k];
	r+=dx*dx;
      }
      r = 1e3*sqrt(r);
      if (r>the_h->r) continue;
      //int64_t type = 0; /* DM */
      //if (po[j].energy>0) type = 1; /* Gas */
      //if (!po[j].energy && po[j].mass < 0.99*PARTICLE_MASS) type = 2; /* Star */
      int64_t type = po[j].type;
      n[type]++;
      m[type]+=po[j].mass;
      check_grp_ids_space(po[j].id+1);
      if (po[j].hid > -1) {
	if (grp_ids[po[j].id]==-1 ||
	    grp_masses[po[j].id] > the_h->m) {
	  grp_ids[po[j].id] = po[j].hid;
	  grp_masses[po[j].id] = the_h->m;
	}
      }
    }
    m[0] /= h0;
    m[1] /= h0;
    m[2] /= h0;

    char *contam = "clean";
    if (m[RTYPE_DM] > 1.1*n[RTYPE_DM]*dm_mass) contam = "contam";
    fprintf(stats, "%-10"PRId64" %-9"PRId64" %-9"PRId64" %-9"PRId64" %-9"PRId64"    %12e  %10f   %12e     %12e    %12e    %11f %11f  %11f   %9f  %9f  %9f  %9f  %9f  %9f  %s  %-10"PRId64"  no\n",
	    the_h->id, n[0]+n[1]+n[2], n[1], n[2], n[0], m[0]+m[1]+m[2], the_h->r/h0, m[1], m[2], m[0], the_h->vmax, the_h->rvmax/h0, the_h->vrms, the_h->pos[0]/h0, the_h->pos[1]/h0, the_h->pos[2]/h0,
	    the_h->pos[3], the_h->pos[4], the_h->pos[5], contam, parents[the_h->id]);
    /*if (the_h->m>1.3e11) {
      printf("%"PRId64" %"PRId64" %"PRId64"\n", n[0], n[1], n[2]);
      }*/
  }
}

void calc_p_start(void) {
  int64_t i;
  for (i=0; i<num_h; i++) { h[i].p_start = -1; h[i].num_p = 0; }
  for (i=0; i<num_p; i++) {
    if (h[p[i].hid].p_start < 0) h[p[i].hid].p_start = i;
    h[p[i].hid].num_p++;
  }
}

int main(int argc, char **argv) {
  int64_t i, j;
  FILE *input;
  char buffer[1024], *end;
  int64_t id, pid;

  if (argc < 5) {
    printf("Usage: %s out.parents tipsy_file.iord highres_dm_mass file1.particles ...\n", argv[0]);
    exit(1);
  }
  
  input = check_fopen(argv[1], "r");
  num_h=0;
  while (fgets(buffer, 1024, input)) {
    if (buffer[0]=='#') continue;
    id = atoll(buffer);
    if (id < 0) continue;
    end = buffer+strlen(buffer);
    while (*end != ' ' && end > buffer) end--;
    pid = atoll(end+1);
    if (id+1 > num_h) {
      num_h = id+1;
      check_realloc_s(parents, sizeof(int64_t), num_h);
    }
    if (pid==id) pid = -1;
    parents[id] = pid;
  }
  fclose(input);
  for (i=0; i<num_h; i++) assert(parents[i] > -2 && parents[i] < num_h);
  input = check_fopen(argv[2], "r"); //Check that it can be opened early
  fclose(input);

  end = argv[1]+strlen(argv[1]);
  while (end>argv[1] && *end!='.') end--;
  if (end>argv[1]) *end='\0';
  snprintf(buffer, 1024, "%s.stat", argv[1]);
  stats=check_fopen(buffer, "w");
  snprintf(buffer, 1024, "%s.grp", argv[1]);
  grp=check_fopen(buffer, "w");
  snprintf(buffer, 1024, "%s.ids", argv[1]);
  ids=check_fopen(buffer, "w");
  fprintf(stats, "  Grp      N_tot     N_gas    N_star    N_dark       Mvir(M_sol)   Rvir(kpc)    GasMass(M_sol)   StarMass(M_sol)   DarkMass(M_sol)       V_max    R@V_max     VelDisp         Xc         Yc         Zc         VXc         VYc         VZc   Contam  Satellite?   False?\n");

  double dm_mass = atof(argv[3]);
  for (i=4; i<argc; i++) {
    num_h = num_p = 0;
    load_full_particles(argv[i], &h, &num_h, &p, &num_p, bounds);
    printf("PARTICLE_MASS: %e\n", PARTICLE_MASS);
    calc_p_start();
    calc_potentials(dm_mass);
  }
  fclose(stats);

  if (strstr(argv[2], ".iord")) {
    input = check_fopen(argv[2], "r");
    int64_t count=0, maxcount=0;
    fgets(buffer, 1024, input);
    maxcount = atoll(buffer);
    if (maxcount > 0) {
      fprintf(grp, "%"PRId64"\n", maxcount);
      fprintf(ids, "%"PRId64"\n", maxcount);
      while (fgets(buffer, 1024, input)) {
	i = atoll(buffer);
	int64_t id = -1;
	if (i>=0 && i<num_grp_ids) id = grp_ids[i];
	fprintf(grp, "%"PRId64"\n", id);
	fprintf(ids, "%"PRId64"\n", i);
	count++;
      }
      if (maxcount != count) {
	fprintf(stderr, "Warning!  Number of IDs actually read (%"PRId64") does not correspond to number of IDs in iord file!  (%"PRId64")\n", count, maxcount);
      }
    } else { //Binary format
      rewind(input);
      struct stat st = {0};
      fstat(fileno(input), &st);
      int64_t n_expected = ((int64_t)st.st_size / (int64_t)4) - (int64_t)1;

      int32_t c=0, mc=0, ii=0, se=0;
      fread(&mc, sizeof(int32_t), 1, input);
      if (mc != n_expected) {
	fprintf(stderr, "Warning: max # of particles (%"PRId32") does not match file length, switching endianness\n", mc);
	swap_endian_4byte((int8_t *)&mc);
	se=1;
      }
      if (mc != n_expected) {
	fprintf(stderr, "Error: max # of particles (%"PRId32") does not match file length (expected: %"PRId64").\n", mc, n_expected);
      }
      if (mc > 0) {
	fprintf(grp, "%"PRId32"\n", mc);
	fprintf(ids, "%"PRId32"\n", mc);
	for (c=0; c<mc; c++) {
	  fread(&ii, sizeof(int32_t), 1, input);
	  if (se) swap_endian_4byte((int8_t *)&ii);
	  int64_t id = -1;
	  if (ii>=0 && ii<num_grp_ids) id = grp_ids[ii];
	  fprintf(grp, "%"PRId64"\n", id);
	  fprintf(ids, "%"PRId32"\n", ii);
	}
      }
    }
    fclose(input);
  } else {
    struct nchilada_header nh;
    char *types[3] = {"gas", "dark", "star"};
    fprintf(grp, "XXXXXXXXXXXXXXXXXX\n");
    fprintf(ids, "XXXXXXXXXXXXXXXXXX\n");
    int64_t count = 0;
    for (j=0; j<3; j++) {
      snprintf(buffer, 1024, "%s/%s/iord", argv[2], types[j]);
      input = check_fopen(buffer, "r");
      int64_t to_read = read_nchilada_header(input, buffer, &nh, 1);
      int64_t id;
      int32_t id32;
      if (nh.code == 5 || nh.code == 6) NCHILADA_ID_BYTES = 4;
      else NCHILADA_ID_BYTES = 8;
      check_fseeko(input, NCHILADA_ID_BYTES*2, SEEK_CUR);
      for (i=0; i<to_read; i++) {
	if (NCHILADA_ID_BYTES == 4) {
	  fread_mswap(&id32, NCHILADA_ID_BYTES, 1, input, NCHILADA_SWAP);
	  id = id32;
	} else {
	  fread_mswap(&id, NCHILADA_ID_BYTES, 1, input, NCHILADA_SWAP);
	}
	fprintf(ids, "%"PRId64"\n", id);
	if (id>=0 && id<num_grp_ids) id = grp_ids[id];
	else id = -1;
	fprintf(grp, "%"PRId64"\n", id);
      }
      count += to_read;
    }
    check_fseeko(grp, 0, SEEK_SET);
    fprintf(grp, "%-18"PRId64"\n", count);
    check_fseeko(ids, 0, SEEK_SET);
    fprintf(ids, "%-18"PRId64"\n", count);
  }
  fclose(grp);
  fclose(ids);
  return(0);
}

