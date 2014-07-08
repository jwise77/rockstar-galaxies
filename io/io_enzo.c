/*
 ENZO INPUT FOR Rockstar
 Modified from Subfind code, courtesy of John Wise
*/
#ifdef ENABLE_HDF5

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>
#include <assert.h>
#include "../universal_constants.h"
#include "../check_syscalls.h"
#include "../config_vars.h"
#include "../rockstar.h"
#include "io_hdf5.h"
#include "io_arepo.h"
#include "io_util.h"

struct EnzoGrid {
  char *fname;
  int64_t level, num_p, id;
};

double EnzoVelocityUnit = 1;
double EnzoMassUnit = 1;
int64_t EnzoNumFiles = 0;

int sort_enzo_grids(const void *a, const void *b) {
  const struct EnzoGrid *c = a;
  const struct EnzoGrid *d = b;
  if (!c->fname) {
    if (!d->fname) return 0;
    return -1;
  }
  if (!d->fname) return 1;
  int res = strcmp(c->fname, d->fname);
  if (res) return res;
  if (c->id < d->id) return -1;
  if (c->id > d->id) return 1;
  fprintf(stderr, "[Warning] Encountered two grids with same id %"PRId64"!\n", c->id);
  return 0;
}

/************************************************************************/
void enzo_load_cosmology (char *filename)
{

  FILE 	*input;
  char 	 buffer[1024];
  float  redshift, initialRedshift;
  int64_t 	 dummy;
  int64_t    TopGrid[3];
  int32_t CycleNumber;
  double Time;
  double Ob = 0.06; //Default for ENZO

  /*********** Get parameters from parameter file ***********/

  input = check_fopen(filename, "r");

  while (fgets(buffer, 1024, input)) {
    sscanf(buffer, "InitialTime        = %lf", &Time);
    sscanf(buffer, "InitialCycleNumber = %d", &CycleNumber);
    sscanf(buffer, "TopGridDimensions   = %"SCNd64" %"SCNd64" %"SCNd64, 
	   TopGrid, TopGrid+1, TopGrid+2);
    sscanf(buffer, "CosmologyHubbleConstantNow = %lf", &h0);
    sscanf(buffer, "CosmologyOmegaMatterNow = %lf", &Om);
    sscanf(buffer, "CosmologyOmegaLambdaNow = %lf", &Ol);
    sscanf(buffer, "CosmologyOmegaBaryonNow = %lf", &Ob);
    sscanf(buffer, "CosmologyComovingBoxSize = %lf", &BOX_SIZE);
    sscanf(buffer, "CosmologyInitialRedshift = %f", &initialRedshift);
    sscanf(buffer, "CosmologyCurrentRedshift = %f", &redshift);
    sscanf(buffer, "#DataCGSConversionFactor[3] = %lg %*s", &EnzoVelocityUnit);
  }  // END line read

  fclose(input);

  /********** Convert to GADGET units **********/

  // rho_crit * omega_M * (comoving_boxsize / rootgrid_reso)^3
  EnzoMassUnit = PARTICLE_MASS = Om * CRITICAL_DENSITY * pow(BOX_SIZE, 3) /
    (TopGrid[0] * TopGrid[1] * TopGrid[2]);

  if (RESCALE_PARTICLE_MASS) EnzoMassUnit *= Om / (Om - Ob + 1e-20);
  AVG_PARTICLE_SPACING = cbrt(PARTICLE_MASS / (Om*CRITICAL_DENSITY));
  SCALE_NOW = 1.0 / (1.0 + redshift);    // Time = a

  /******* Get number of files and particles from .hierarchy file *******/
  snprintf(buffer, 1024, "%s.hierarchy", filename);
  EnzoNumFiles = 0;
  
  input = check_fopen(buffer, "r");
  //printf("enzoFindFiles: parsing hierarchy file ...\n"); 

  while(fgets(buffer, 1024, input)) {
    if (sscanf(buffer, "NumberOfParticles = %"PRId64"\n", &dummy) == 1)
      EnzoNumFiles++;
  }  // END file line read
  fclose(input);
  //fprintf(stderr, "Grids = %"PRId64"\n", EnzoNumFiles);
}

/************************************************************************/

void load_particles_enzo(char *filename, struct particle **p, int64_t *num_p) {

  char buffer[1024];
  int64_t i, dim, dummy, grid;
  struct EnzoGrid *grids = NULL;
  int64_t startIndex, endIndex;
  float GridLeftEdge[3], GridRightEdge[3];
  float	cellWidth, rootCellWidth=1;

  int64_t length = strlen(filename);
  int64_t block = 0;

  if (NUM_BLOCKS > 1) {
    length--;
    while (length && filename[length] != '.') length--;
    assert(filename[length]=='.');
    filename[length] = 0;
    block = atol(filename+length+1);
  }

  enzo_load_cosmology(filename);
  check_realloc_s(grids, sizeof(struct EnzoGrid), EnzoNumFiles);
  memset(grids, 0, sizeof(struct EnzoGrid)*EnzoNumFiles);

  /******** Get particle counts and levels in each grid ********/
  snprintf(buffer, 1024, "%s.hierarchy", filename);
  FILE *hf = check_fopen(buffer, "r");
  //printf("load_particles_enzo: parsing hierarchy file ...\n");

  grid = 0;
  while (fgets(buffer, 1024, hf)) {
    if (sscanf(buffer, "Grid = %"SCNd64, &dummy) == 1) {
      grid = dummy-1;
      grids[grid].id = dummy;
    }
    sscanf(buffer, "GridStartIndex = %"SCNd64, &startIndex);
    sscanf(buffer, "GridEndIndex = %"SCNd64, &endIndex);
    sscanf(buffer, "GridLeftEdge = %f %f %f", GridLeftEdge+0, GridLeftEdge+1, 
	   GridLeftEdge+2);
    sscanf(buffer, "GridRightEdge = %f %f %f", GridRightEdge+0, GridRightEdge+1,
	   GridRightEdge+2);
    if (sscanf(buffer, "NumberOfParticles = %"SCNd64, &dummy) == 1) {
      // Calculate and store grid level
      cellWidth = (GridRightEdge[0] - GridLeftEdge[0]) / (endIndex - startIndex + 1);
      if (grid == 0) rootCellWidth = cellWidth;
      grids[grid].level = log2(rootCellWidth / cellWidth) + 0.5;
      grids[grid].num_p = dummy;
    }

    if (!strncmp(buffer, "ParticleFileName = ", 19)) {
      grids[grid].fname = strdup(buffer + 19);
      char *end = grids[grid].fname+strlen(grids[grid].fname);
      if (end > grids[grid].fname && end[-1]=='\n')
	end[-1] = 0;
    }
  }
  fclose(hf);

  qsort(grids, EnzoNumFiles, sizeof(struct EnzoGrid), sort_enzo_grids);

  // Particle HDF labels
  char *ParticlePositionLabel[] = 
    {"particle_position_x", "particle_position_y", "particle_position_z"};
  char *ParticleVelocityLabel[] =
    {"particle_velocity_x", "particle_velocity_y", "particle_velocity_z"};
  char *ParticleMassLabel = "particle_mass";
  char *ParticleIDLabel = "particle_index";
  //char *ParticleTypeLabel = "particle_type";

  int64_t total_p = 0, p_start, to_read, total_p2 = 0;
  for (i=0; i<EnzoNumFiles; i++) total_p += grids[i].num_p;
  particle_range(total_p, block, NUM_BLOCKS, &p_start, &to_read);
  int64_t b_start = -1, b_end = -1;
  for (i=0; i<EnzoNumFiles; i++) {
    if (b_start < 0 && total_p2 >= p_start) b_start = i;
    total_p2 += grids[i].num_p;
    if (b_end < 0 && total_p2 >= p_start+to_read) b_end = i+1;
  }
  assert(b_start > -1 && b_end > -1);  

  char *lastfile = 0;
  hid_t h5_file = 0;
  int64_t total_alloced = *num_p;
  int64_t total_read = 0;
  for(i=b_start; i<b_end; i++) {
    if (b_start == b_end) break;
    if (!grids[i].fname) continue;
    snprintf(buffer, 1024, "Grid%8.8"PRId64, grids[i].id);
    if (!lastfile || strcmp(lastfile, grids[i].fname)!=0) {
      if (lastfile) H5Fclose(h5_file);
      lastfile = grids[i].fname;
      h5_file = check_H5Fopen(lastfile, H5F_ACC_RDONLY);
      //printf("Opened %s!\n", grids[i].fname);
    }
      
    check_realloc_smart(*p, sizeof(struct particle), total_alloced, 
			(*num_p) + grids[i].num_p);
    memset(*p + (*num_p), 0, sizeof(struct particle)*grids[i].num_p);

    for (dim = 0; dim < 3; dim++) {
      arepo_read_dataset(h5_file, lastfile, buffer, ParticlePositionLabel[dim],
			 *p + (*num_p), grids[i].num_p, 
			 (char *)&(p[0][0].pos[dim])-(char*)(p[0]),
			 1, H5T_NATIVE_FLOAT);
      arepo_read_dataset(h5_file, lastfile, buffer, ParticleVelocityLabel[dim],
			 *p + (*num_p), grids[i].num_p,
			 (char *)&(p[0][0].pos[dim+3])-(char*)(p[0]),
			 1, H5T_NATIVE_FLOAT);
    }

    arepo_read_dataset(h5_file, lastfile, buffer, ParticleMassLabel,
		       *p + (*num_p), grids[i].num_p, 
		       (char *)&(p[0][0].mass)-(char*)(p[0]),
		       1, H5T_NATIVE_FLOAT);
    arepo_read_dataset(h5_file, lastfile, buffer, ParticleIDLabel,
		       *p + (*num_p), grids[i].num_p,
		       (char *)&(p[0][0].id)-(char*)(p[0]),
		       1, H5T_NATIVE_LLONG);

    for(int64_t n=(*num_p); n<(*num_p)+grids[i].num_p; n++) {
      for (dim = 0; dim < 3; dim++) {
	p[0][n].pos[dim] *= BOX_SIZE;
	p[0][n].pos[dim+3] *= EnzoVelocityUnit / 1.0e5; //CGS to km/s
      }

      p[0][n].mass *= EnzoMassUnit / pow(8.0, grids[i].level);
      p[0][n].type = RTYPE_DM;
    }
    (*num_p) += grids[i].num_p;
    total_read += grids[i].num_p;
  }

  if (lastfile) H5Fclose(h5_file);

  for (i=0; i<EnzoNumFiles; i++) 
    if (grids[i].fname) free(grids[i].fname);
  free(grids);
  //printf("Particles read: %"PRId64"\n", total_read);

  if (NUM_BLOCKS > 1) filename[length] = '.';
  snprintf(buffer, 1024, "particles.%"PRId64".ascii", block);
  FILE *output = check_fopen(buffer, "w");
  fprintf(output, "#X Y Z VX VY VZ Mass Energy ID Type\n");
  for (i=0; i<*num_p; i++) {
    fprintf(output, "%f %f %f %f %f %f %f %f %"PRId64" %d\n",
	    p[0][i].pos[0], p[0][i].pos[1], p[0][i].pos[2], 
	    p[0][i].pos[3], p[0][i].pos[4], p[0][i].pos[5], 
	    p[0][i].mass, p[0][i].energy, p[0][i].id, p[0][i].type); 
  }
  fclose(output);
}


#endif /*ENABLE_HDF5*/
