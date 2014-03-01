/* From Matt Turk, Copyright (c) 2011.
   Explicit permission granted to distribute with Rockstar under the GPLv3.
*/
#ifndef _IO_GENERIC_H
#define _IO_GENERIC_H

#include <stdint.h>
#include "../particle.h"
#include "../potential.h"
#include "../halo.h"

typedef void (*LPG) (char *filename, struct particle **p, int64_t *num_p);
typedef void (*AHG) (struct halo *h, struct potential *p, int64_t num_p);

extern LPG load_particles_generic;
extern AHG analyze_halo_generic;

void set_load_particles_generic(LPG func, AHG afunc);

#endif /* _IO_GENERIC_H */
