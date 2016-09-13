#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <rpc/types.h>
#include <rpc/xdr.h>
#include "../check_syscalls.h"
#include "../universal_constants.h"
#include "../config_vars.h"
#include "../config.h"
#include "../particle.h"
#include "io_tipsy.h"

enum tipsy_DataTypeCode {
	int8 = 1,
	uint8,
	int16,
	uint16,
	int32,
	uint32,
	int64,
	uint64,
	float32,
	float64
};


void load_particles_tipsy(char *filename, struct particle **p, int64_t *num_p) {
  FILE *input;
  struct tipsy_dump header;
  struct tipsy_dark_particle dark;
  struct tipsy_star_particle star;
  struct tipsy_gas_particle gas;

  int i, j;
  int xdrfmt=1,haveiords=0;
  int *iords=NULL;
  XDR xdrs;

  input = check_fopen(filename, "r");
  xdrstdio_create(&xdrs, input, XDR_DECODE);
  tipsy_xdr_header(&xdrs, &header);
  if (header.ndim != 3) {
      xdr_destroy(&xdrs);
      fseek(input, 0L, SEEK_SET);
      check_fread((char *)&header,sizeof(header),1,input);
      assert(header.ndim == 3);
      xdrfmt = 0;
      }
  SCALE_NOW = header.time;
  printf("Gas particles: %d\nDark Matter Particles: %d\nStar Particles: %d\n",
	 header.nsph, header.ndark, header.nstar);
  check_realloc_s((*p), ((*num_p) + header.ndark+header.nsph+header.nstar),
		  sizeof(struct particle));

  haveiords = load_ids_tipsy(filename,header,&iords);

  for(i = 0; i < header.nsph; i++) {
    if (xdrfmt) assert(tipsy_xdr_gas(&xdrs, &gas) > 0);
    else fread((char *)&gas,sizeof(struct tipsy_gas_particle), 1, input) ;
    for (j=0; j<3; j++) {
      if (haveiords) (*p)[i+(*num_p)].id = iords[i];
      else (*p)[i+(*num_p)].id = i;
      (*p)[i+(*num_p)].pos[j] = (gas.pos[j] + 0.5) * TIPSY_LENGTH_CONVERSION;
      (*p)[i+(*num_p)].pos[j+3] = gas.vel[j] * TIPSY_VELOCITY_CONVERSION *SCALE_NOW;
    }
    (*p)[i+(*num_p)].type = RTYPE_GAS;
    (*p)[i+(*num_p)].mass = gas.mass * TIPSY_MASS_CONVERSION;
    /* Specific energy! */
    (*p)[i+(*num_p)].energy = gas.temp * ENERGY_PER_KELVIN_PER_UNIT_MASS;
  }
  (*num_p) += header.nsph;  

  for(i = 0; i < header.ndark; i++) {
    int ip = i+header.nsph;
    if (xdrfmt) assert(tipsy_xdr_dark(&xdrs, &dark) > 0);
    else check_fread((char *)&dark,sizeof(struct tipsy_dark_particle), 1, input) ;
    for (j=0; j<3; j++) {
      if (haveiords) (*p)[i+(*num_p)].id = iords[ip];
      else (*p)[i+(*num_p)].id = ip;
      (*p)[i+(*num_p)].pos[j] = (dark.pos[j] + 0.5) * TIPSY_LENGTH_CONVERSION;
      (*p)[i+(*num_p)].pos[j+3] = dark.vel[j] * TIPSY_VELOCITY_CONVERSION*SCALE_NOW;
    }
    (*p)[i+(*num_p)].type = RTYPE_DM;
    (*p)[i+(*num_p)].mass = dark.mass * TIPSY_MASS_CONVERSION;
    (*p)[i+(*num_p)].energy = 0;
  }

  //printf("Read %d dark matter particles.\n",header.ndark);
  //TOTAL_PARTICLES += header.ndark;
  (*num_p) += header.ndark;

  for(i = 0; i < header.nstar; i++) {
    int ip = i+header.nsph+header.ndark;
    if (xdrfmt) assert(tipsy_xdr_star(&xdrs, &star) > 0);
    else fread((char *)&star,sizeof(struct tipsy_star_particle), 1, input);
    for (j=0; j<3; j++) {
      if (haveiords) (*p)[i+(*num_p)].id = iords[ip];
      else (*p)[i+(*num_p)].id = ip;
      (*p)[i+(*num_p)].pos[j] = (star.pos[j] + 0.5) * TIPSY_LENGTH_CONVERSION;
      (*p)[i+(*num_p)].pos[j+3] = star.vel[j] * TIPSY_VELOCITY_CONVERSION * SCALE_NOW;
    }
    (*p)[i+(*num_p)].type = RTYPE_STAR;
    (*p)[i+(*num_p)].mass = star.mass * TIPSY_MASS_CONVERSION;
    (*p)[i+(*num_p)].energy = 0;
  }
  (*num_p) += header.nstar;

  if (haveiords) free(iords);
  if (xdrfmt) xdr_destroy(&xdrs);
  fclose(input);
}

/* open iord file for ids */
int load_ids_tipsy(char *filename, struct tipsy_dump header, int **iords) {
  FILE *iordf;
  XDR xdrs;
  char iofilename[256];
  int i, nbodies=0, count=0, bStandard = 0, bASCII=0;
  sprintf(iofilename,"%s.iord",filename);
  iordf = fopen(iofilename, "r");
  if (iordf != NULL) {
      count=fscanf(iordf, "%d%*[, \t]%*d%*[, \t]%*d",&nbodies) ;
      if ( (count == EOF) || (count==0) ){
	  /* try binary instead */
	  rewind(iordf);
	  count=fread(&nbodies, sizeof(int), 1, iordf) ;

	  if ( (count == EOF) || (count==0) ){
	      printf("<%s format is wrong>\n",iofilename);
	      fclose(iordf);
	      return 0;
	  } else if (nbodies != header.nbodies) {
	      fseek(iordf,0,SEEK_SET);
	      xdrstdio_create(&xdrs,iordf,XDR_DECODE);
	      xdr_int(&xdrs,&nbodies);
	      if (nbodies != header.nbodies) {
		printf("<%s doesn't appear standard or binary or nbodies (%d) != expected (%d).>\n",iofilename, nbodies, header.nbodies);
		  xdr_destroy(&xdrs);
		  fclose(iordf);
		  return 0;
		  } else bStandard = 1;
	      }
	  } else bASCII = 1;

      /* allocate iords array */
      if(*iords != NULL) free(*iords);
      *iords = (int *)malloc(nbodies*sizeof(**iords));
      if(*iords == NULL) {
	  printf("no room for iords\n");
	  fclose(iordf);
	  return 0;
	  }

      for(i = 0, count = 0; i < nbodies; i++) {
	int idummy, check=0;
	  if(count >= nbodies) {
	      printf("<%s format is wrong>\n",iofilename);
	      free(*iords) ;
	      *iords = NULL;
	      return 0;
	      }
	  if (bStandard) check = xdr_int(&xdrs,&idummy);
	  else if (bASCII) check = fscanf(iordf, "%d", &idummy);
	  else check = check_fread(&idummy, sizeof(int), 1, iordf) ;
	  (*iords)[count] = idummy;

	  if(check == EOF) {
	      printf("<%s format is wrong>\n",iofilename);
	      free(*iords) ;
	      *iords = NULL;
	      break;
	      }
	  count++;
	  }
      fclose(iordf);
      }
  if (iordf == NULL || count != nbodies) {
    printf("WARNING: Did not read iords file\n");
    return 0;
  }
  //printf("Read %d iords.\n",count);
  return 1;
}

int tipsy_xdr_header(XDR *pxdrs,struct tipsy_dump *ph) {
    int pad = 0;
    
    if (!xdr_double(pxdrs,&ph->time)) return 0;
    if (!xdr_int(pxdrs,&ph->nbodies)) return 0;
    if (!xdr_int(pxdrs,&ph->ndim)) return 0;
    if (!xdr_int(pxdrs,&ph->nsph)) return 0;
    if (!xdr_int(pxdrs,&ph->ndark)) return 0;
    if (!xdr_int(pxdrs,&ph->nstar)) return 0;
    if (!xdr_int(pxdrs,&pad)) return 0;
    return 1;
    }

int tipsy_xdr_gas(XDR *pxdrs,struct tipsy_gas_particle *ph) { 
    int i;
    if (!xdr_float(pxdrs,&ph->mass)) return 0;
    for(i=0; i<TIPSY_MAXDIM; i++) if (!xdr_float(pxdrs,&ph->pos[i])) return 0;
    for(i=0; i<TIPSY_MAXDIM; i++) if (!xdr_float(pxdrs,&ph->vel[i])) return 0;
    if (!xdr_float(pxdrs,&ph->rho)) return 0;
    if (!xdr_float(pxdrs,&ph->temp)) return 0;
    if (!xdr_float(pxdrs,&ph->hsmooth)) return 0;
    if (!xdr_float(pxdrs,&ph->metals)) return 0;
    if (!xdr_float(pxdrs,&ph->phi)) return 0;
    return 1;
    }

int tipsy_xdr_dark(XDR *pxdrs,struct tipsy_dark_particle *ph) {
    int i;
    if (!xdr_float(pxdrs,&ph->mass)) return 0;
    for(i=0; i<TIPSY_MAXDIM; i++) if (!xdr_float(pxdrs,&ph->pos[i])) return 0;
    for(i=0; i<TIPSY_MAXDIM; i++) if (!xdr_float(pxdrs,&ph->vel[i])) return 0;
    if (!xdr_float(pxdrs,&ph->eps)) return 0;
    if (!xdr_float(pxdrs,&ph->phi)) return 0;
    return 1;
    }

int tipsy_xdr_star(XDR *pxdrs,struct tipsy_star_particle *ph) {
    int i;
    if (!xdr_float(pxdrs,&ph->mass)) return 0;
    for(i=0; i<TIPSY_MAXDIM; i++) if (!xdr_float(pxdrs,&ph->pos[i])) return 0;
    for(i=0; i<TIPSY_MAXDIM; i++) if (!xdr_float(pxdrs,&ph->vel[i])) return 0;
    if (!xdr_float(pxdrs,&ph->metals)) return 0;
    if (!xdr_float(pxdrs,&ph->tform)) return 0;
    if (!xdr_float(pxdrs,&ph->eps)) return 0;
    if (!xdr_float(pxdrs,&ph->phi)) return 0;
    return 1;
    }

