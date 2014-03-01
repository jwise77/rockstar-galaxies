/* From Matt Turk, Copyright (c) 2011.
   Explicit permission granted to distribute with Rockstar under the GPLv3.
*/
#include <stdlib.h>
#include <string.h>
#include "io_generic.h"

LPG load_particles_generic = NULL;
AHG analyze_halo_generic = NULL;

void set_load_particles_generic(LPG func, AHG afunc)
{
  load_particles_generic = func;
  analyze_halo_generic = afunc;
}
