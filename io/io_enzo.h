#ifndef _IO_ENZO_H_
#define _IO_ENZO_H_
#ifdef ENABLE_HDF5
#include <stdint.h>
#include "../particle.h"

void load_particles_enzo(char *filename, struct particle **p, int64_t *num_p);

#endif /* ENABLE_HDF5 */
#endif /* _IO_ENZO_H_ */
