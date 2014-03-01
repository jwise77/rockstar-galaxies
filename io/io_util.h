#ifndef _IO_UTIL_H_
#define _IO_UTIL_H_

#include <stdlib.h>
#include <stdint.h>

void swap_endian_4byte(int8_t *e);
void swap_endian_8byte(int8_t *e);
void swap_4byte_to_8byte(int32_t *e);
size_t fread_swap(void *ptr, size_t size, size_t nitems, FILE *stream);
size_t fread_swap8(void *ptr, size_t size, size_t nitems, FILE *stream);
size_t fread_fortran(void *ptr, size_t size, size_t nitems, FILE *stream, int swap);
size_t fwrite_fortran(void *ptr, size_t size, size_t nitems, FILE *stream);
void skip_fortran(FILE *stream, int swap);
void particle_range(int64_t total_p, int64_t block, int64_t num_blocks,
		    int64_t *p_start, int64_t *to_read);
size_t fread_mswap(void *ptr, size_t size, size_t nitems, FILE *stream, int swap);

#endif /* _IO_UTIL_H_ */
