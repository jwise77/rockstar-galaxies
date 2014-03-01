#ifndef _IO_ART_H_
#define _IO_ART_H_

#define ART_TEXT_SIZE 45
#define EXTRA_VARS 100
#include <stdint.h>
#include "../particle.h"

struct PMcrd_header {
  char header[ART_TEXT_SIZE];
  float AEXPN,AEXP0,AMPLT,ASTEP;
  int32_t ISTEP;
  float PARTW,TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0;
  int32_t NROWC,NGRIDC,Nspecies,Nseed;
  float Om0,Oml0,hubble,Wp5,Ocurv,extras[EXTRA_VARS-1], Box;
} __attribute__ ((packed));

struct art_header1 {
  float AEXPN, ASTEP;
  int32_t ISTEP, NROWC, NGRIDC, Nspecies, Nseed;
  float Om0, Oml0, hubble, Box;
} __attribute__ ((packed));

struct art_header1a {
  float AEXPN, ASTEP;
  int32_t ISTEP, NROWC, NGRIDC, Nspecies, Nseed;
  float Om0, Oml0, hubble, Box, MassOne;
} __attribute__ ((packed));

struct art_header2 {
  int32_t k, Nx, Ny, Nz;
  float dR;
} __attribute__ ((packed));

struct art_header2a {
  int32_t node, Nx, Ny, Nz;
  float dBuffer;
  int32_t nBuffer;
} __attribute__ ((packed));

struct art_header3 {
  float xL,xR,yL,yR,zL,zR;
} __attribute__ ((packed));

struct art_particle {
  float pos[6];
  uint64_t id;
} __attribute__ ((packed));

void load_particles_art(char *filename, struct particle **p, int64_t *num_p);

#endif /* _IO_ART_H_ */
