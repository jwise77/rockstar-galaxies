#ifndef PARTICLE_H
#define PARTICLE_H
#include <stdint.h>

#define RTYPE_DM   0
#define RTYPE_GAS  1
#define RTYPE_STAR 2
#define RTYPE_BH   3
#define NUM_RTYPES 4

struct particle {
  int64_t id;
  float pos[6];
  float mass, energy; /*Energy per unit mass*/
  float softening; /* Per-particle softening, not currently used. */
  float metallicity; /* Not currently used. */
  int32_t type;
};

#endif /*PARTICLE_H*/
