#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>
#include <string.h>
#include <strings.h>
#include <sys/resource.h>
#include "config_vars.h"
#include "check_syscalls.h"
#include "io/read_config.h"
#include "universal_constants.h"

void setup_config(void) {
  if (!NUM_READERS) NUM_READERS = NUM_BLOCKS;

  if (!PARTICLE_MASS)
    PARTICLE_MASS = CRITICAL_DENSITY * BOX_SIZE * BOX_SIZE * BOX_SIZE
      * Om / TOTAL_PARTICLES;
  if (!AVG_PARTICLE_SPACING)
    AVG_PARTICLE_SPACING = cbrt(PARTICLE_MASS / (Om*CRITICAL_DENSITY));

  if (LIGHTCONE || !PARALLEL_IO) {
    PERIODIC = 0;
    TEMPORAL_HALO_FINDING = 0;
  }

  if (IGNORE_PARTICLE_IDS)
    TEMPORAL_HALO_FINDING = 0;

  if (!FORCE_RES_PHYS_MAX)
    FORCE_RES_PHYS_MAX = FORCE_RES;

  struct rlimit rlp;
  getrlimit(RLIMIT_NOFILE, &rlp);
  rlp.rlim_cur = rlp.rlim_max;
  setrlimit(RLIMIT_NOFILE, &rlp);
  getrlimit(RLIMIT_CORE, &rlp);
  rlp.rlim_cur = rlp.rlim_max;
  setrlimit(RLIMIT_CORE, &rlp);
  if (NUM_WRITERS < FORK_PROCESSORS_PER_MACHINE)
    NUM_WRITERS = FORK_PROCESSORS_PER_MACHINE;

  if (STARTING_SNAP >= NUM_SNAPS) {
    fprintf(stderr, "[Warning] No work will be done unless NUM_SNAPS > STARTING_SNAP in config file!\n");
  }
}

void do_config(char *filename) {
  struct configfile c = {0};
  if (filename && strlen(filename))
    load_config(&c, filename);

#define string(a,b) a = config_to_string(&c,#a,b)
#define real(a,b) a = config_to_real(&c,#a,b)
#define real3(a,b) config_to_real3(&c,#a,a,b)
#define integer(a,b) a = config_to_real(&c,#a,b)
#include "config.template.h"
#undef string
#undef real
#undef real3
#undef integer
  syntax_check(&c, "[Warning]");
  setup_config();
  free_config(c);
  if (filename && strlen(filename)) {
    free(ROCKSTAR_CONFIG_FILENAME);
    ROCKSTAR_CONFIG_FILENAME = strdup(filename);
  }
}

void output_config(char *filename) {
  char buffer[1024];
  FILE *output;
  if (!filename) filename = "rockstar.cfg";
  snprintf(buffer, 1024, "%s/%s", OUTBASE, filename);
  output = check_fopen(buffer, "w");

#define string(a,b) fprintf(output, "%s = \"%s\"\n", #a, a);
#define real(a,b) fprintf(output, "%s = %g\n", #a, a);
#define real3(a,b) fprintf(output, "%s = (%g, %g, %g)\n", #a, a[0], a[1], a[2]);
#define integer(a,b) fprintf(output, "%s = %"PRId64"\n", #a, a);
#include "config.template.h"

  fclose(output);
}
