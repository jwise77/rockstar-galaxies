#ifndef _IO_ASCII_H
#define _IO_ASCII_H

#include <stdint.h>
#include "../particle.h"

void gzip_file(char *filename);
void load_particles(char *filename, struct particle **p, int64_t *num_p);

#endif /* _IO_ASCII_H */
