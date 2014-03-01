#ifndef _BIT_ARRAY_H_
#define _BIT_ARRAY_H_
#include <inttypes.h>
#include <string.h>
#include "check_syscalls.h"

#define CMAX ((1<<8)-1)

#define BIT_ALLOC(n) (char *)check_realloc(NULL, ((n)/(int64_t)8)+1, "bit array");

//Note that X must be a character array!
#define BIT_ALL_CLEAR(x,n) memset(x, 0, ((n)/(int64_t)8)+1)

//Set bit Y in bitarray X to 1
#define BIT_SET(x,y) (x[(uint64_t)(y)>>((int64_t)3)] |= (1<<((y)&7)))

//Set bit Y in bitarray X to 0
#define BIT_CLR(x,y) (x[(uint64_t)(y)>>((int64_t)3)] &= (CMAX-(1<<((y)&7))))

//Returns nonzero value if bit Y in bitarray X is 1
#define BIT_TST(x,y) (x[(uint64_t)(y)>>((int64_t)3)] & ((1<<((y)&7))))

#endif /* _BIT_ARRAY_H_ */
