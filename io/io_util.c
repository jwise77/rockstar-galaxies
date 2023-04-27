#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <assert.h>
#include <limits.h>
#include "io_util.h"
#include "../check_syscalls.h"

void swap_endian_4byte(int8_t *e) {
  int8_t a;
#define SWAP(x,y) a=e[x]; e[x]=e[y]; e[y]=a;
  SWAP(0,3);
  SWAP(1,2);
}

void swap_endian_8byte(int8_t *e) {
  int8_t a;
  SWAP(0,7);
  SWAP(1,6);
  SWAP(2,5);
  SWAP(3,4);
}

void swap_4byte_to_8byte(int32_t *e) {
  int32_t a;
  SWAP(0,1);
}
#undef SWAP


size_t fread_swap(void *ptr, size_t size, size_t nitems, FILE *stream) {
  size_t i, n;
  int8_t *e = ptr;
  n = check_fread(ptr, size, nitems, stream);
  if (!(size % 4))
    for (i=0; i<(n*size); i+=4) swap_endian_4byte(e+i);
  return n;
}

size_t fread_swap8(void *ptr, size_t size, size_t nitems, FILE *stream) {
  size_t i, n;
  int8_t *e = ptr;
  n = check_fread(ptr, size, nitems, stream);
  assert(!((n*size) % 8));
  for (i=0; i<(n*size); i+=8) swap_endian_8byte(e+i);
  return n;
}

size_t fread_mswap(void *ptr, size_t size, size_t nitems, FILE *stream, int swap) {
  if (swap) {
    if (size==8) return fread_swap8(ptr, size, nitems, stream);
    else return fread_swap(ptr, size, nitems, stream);
  }
  return check_fread(ptr, size, nitems, stream);
}

size_t _fread_fortran(void *ptr, size_t size, size_t nitems, FILE *stream, int swap) {
  uint32_t head, tail;
  size_t n;
  if (swap) fread_swap(&head, 4, 1, stream);
  else check_fread(&head, 4, 1, stream);
  assert((size*nitems) >= head);
  if (swap) n = fread_swap(ptr, size, nitems, stream);
  else n = check_fread(ptr, size, nitems, stream);
  if (swap) fread_swap(&tail, 4, 1, stream);
  else check_fread(&tail, 4, 1, stream);
  if (head != tail) {
    printf("Header length: %d; tail length: %d\n", head, tail);
    assert(head == tail);
  }
  return n;
}

size_t fread_fortran(void *ptr, size_t size, size_t nitems, FILE *stream, int swap) {
  int64_t nread = 0;
  while (nread < nitems)
    nread += _fread_fortran(((char *)ptr) + (size*nread), size, nitems - nread, 
			    stream, swap);
  return ((size_t) nread);
}

void skip_fortran(FILE *stream, int swap) {
  uint32_t head, tail;
  if (swap) fread_swap(&head, 4, 1, stream);
  else check_fread(&head, 4, 1, stream);
  check_fseeko(stream, head, SEEK_CUR);
  if (swap) fread_swap(&tail, 4, 1, stream);
  else check_fread(&tail, 4, 1, stream);
  assert(head == tail);
}

size_t _fwrite_fortran(void *ptr, size_t size, size_t nitems, FILE *stream) {
  uint32_t head;
  assert(size*nitems <= (uint64_t)(UINT32_MAX));
  head = size*nitems;
  check_fwrite(&head, sizeof(uint32_t), 1, stream);
  check_fwrite(ptr, size, nitems, stream);
  check_fwrite(&head, sizeof(uint32_t), 1, stream);
  return nitems;
}

size_t fwrite_fortran(void *ptr, size_t size, size_t nitems, FILE *stream) {
  int64_t nwritten = 0, n_to_write;
  assert(size > 0);
  while (nwritten < nitems) {
    n_to_write = (((uint64_t)UINT32_MAX) / size);
    if (n_to_write > (nitems - nwritten)) n_to_write = nitems - nwritten;
    _fwrite_fortran(((char *)ptr) + (size*nwritten), size, n_to_write, stream);
    nwritten += n_to_write;
  }
  return ((size_t) nwritten);
}

void particle_range(int64_t total_p, int64_t block, int64_t num_blocks,
		    int64_t *p_start, int64_t *to_read) {
  int64_t particles_per_block = (total_p+num_blocks-1) / num_blocks;
  *p_start = particles_per_block * block;
  *to_read = particles_per_block;
  if (*p_start > total_p) *p_start = total_p;
  if (*p_start + *to_read > total_p)
    *to_read = (total_p - *p_start);
}
