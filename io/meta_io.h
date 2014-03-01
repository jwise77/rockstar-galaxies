#ifndef _META_IO_H_
#define _META_IO_H_
#include <stdint.h>
#include "../halo.h"
#include "io_internal.h"

extern char **snapnames;
extern char **blocknames;

int _should_print(struct halo *h, float *bounds);
void read_input_names(char *filename, char ***stringnames, int64_t *num_names);
void get_input_filename(char *buffer, int maxlen, int64_t snap, int64_t block);
void get_output_filename(char *buffer, int maxlen, int64_t snap, 
			 int64_t chunk, char *type);
void read_particles(char *filename);
int64_t count_halos_to_print(float *bounds);
void delete_binary(int64_t snap, int64_t chunk);
void output_halos(int64_t id_offset, int64_t snap, 
		  int64_t chunk, float *bounds);

char *gen_merger_catalog(int64_t snap, int64_t chunk, struct halo *halos, int64_t num_halos, int64_t *cat_length, int64_t *header_length);
void output_merger_catalog(int64_t snap, int64_t chunk, int64_t location, int64_t length, char *cat);

#endif /* _META_IO_H_ */
