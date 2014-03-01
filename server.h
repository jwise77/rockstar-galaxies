#ifndef _SERVER_H_
#define _SERVER_H_

#define READER_TYPE 0
#define WRITER_TYPE 1
#define DATA_PARTICLES 0
#define DATA_BPARTICLES 1
#define DATA_HALOS 2
#define DATA_HPARTICLES 3


#include <stdint.h>

struct client_info {
  int type;
  int64_t cs;
  char *address, *serv_port;
  int port;
  float bounds[6];
  float halo_bounds[6];
  float bounds_prevsnap[6];
  int64_t num_halos;
  int64_t status;
  int64_t workers;
  int64_t head_length;
  int64_t cat_length;
};

#define PROJECTION_SIZE 10000

struct projection_request {
  int64_t dir, id;
  float bounds[6];
};

struct projection {
  int64_t dir, id;
  float bounds[6];
  int64_t data[PROJECTION_SIZE];
};


int server(void);
void check_num_writers(void);

#define timed_output(...) { print_time(); fprintf(stderr, __VA_ARGS__); }

#endif /* _SERVER_H_ */
