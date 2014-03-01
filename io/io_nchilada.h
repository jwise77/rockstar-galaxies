#ifndef _IO_NCHILADA_H_
#define _IO_NCHILADA_H_
#include <inttypes.h>

#define NCHILADA_MAGIC ((int32_t)1062053)

struct nchilada_header {
 int32_t    magic;
 double time;
 int32_t    iHighWord;
 int32_t    nbodies;
 int32_t    ndim;
 int32_t    code;
};

int64_t read_nchilada_header(FILE *input, char *filename, struct nchilada_header *nh, int64_t dims);
void load_particles_nchilada(char *filename, struct particle **p, int64_t *num_p);

#endif /* _IO_NCHILADA_H_ */

