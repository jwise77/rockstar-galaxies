#ifndef _STRINGPARSE_H_
#define _STRINGPARSE_H_
#include <inttypes.h>

enum parsetype {
  PARSE_FLOAT32 = 0,
  PARSE_INT32 = 1,
  PARSE_FLOAT64 = 2,
  PARSE_INT64 = 3,
  PARSE_STRING = 4,
  PARSE_SKIP = 5,
};

#define short_parsetype parsetype

#define SHORT_PARSETYPE \
enum parsetype F; F = PARSE_FLOAT32;		\
enum parsetype D; D = PARSE_INT32;		\
enum parsetype F64; F64 = PARSE_FLOAT64;	\
enum parsetype LF; LF = PARSE_FLOAT64;		\
enum parsetype LD; LD = PARSE_INT64;		\
enum parsetype D64; D64 = PARSE_INT64;		\
enum parsetype S; S = PARSE_STRING;		\
enum parsetype K; K = PARSE_SKIP;		

struct parse_format {
  int64_t index;
  enum parsetype type;
  void *data;
};

int64_t stringparse(char *buffer, void **data, enum parsetype *types, int64_t max_n);
int64_t stringparse_format(char *buffer, struct parse_format *pf, int64_t max_n);

#endif
