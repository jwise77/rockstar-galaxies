#ifndef _RSOCKET_H_
#define _RSOCKET_H_
#include <inttypes.h>

#define RSOCKET_READ 1
#define RSOCKET_WRITE 2
#define RSOCKET_ERROR 4
#define RSOCKET_VERIFY (uint64_t)0xfadedaceabeec0c0

struct rpacket_header {
  int64_t seq;
  int64_t length;
  int64_t flags;
};

struct rsocket {
  int32_t fd;
  int64_t id;
  int64_t server_id;
  uint64_t magic;
  int64_t sseq, rseq;
  int64_t flags;
  void *last_data;
  char *address;
  char *port;
};

int64_t connect_to_addr(char *host, char *port);
int64_t listen_at_addr(char *host, char *port);
void close_rsocket(int64_t s);
int rsocket_fd(int64_t s);
int64_t rsocket_from_fd(int fd);
int64_t accept_connection(int64_t s, char **address, int *port);
int64_t send_to_socket(int64_t s, void *data, int64_t length);
int64_t send_msg(int64_t s, void *data, int64_t length);
int64_t recv_from_socket(int64_t s, void *data, int64_t length);
void *recv_and_alloc(int64_t s, void *data, int64_t length);
void *recv_msg(int64_t s, void *data, int64_t *length, int64_t offset);
void *recv_msg_nolength(int s, void *data);
int64_t send_to_socket_noconfirm(int64_t s, void *data, int64_t length);
int64_t send_to_socket_delayconfirm(int64_t s, void *data, int64_t length);
int64_t send_to_socket_confirm(int64_t s, void *data, int64_t length);

void tag_rsocket(int64_t s);
void clear_rsocket_tags(void);
int64_t check_rsocket_tag(int64_t s);
int64_t select_rsocket(int select_type, double timeout);

//Internal
int64_t rsocket_accept_connection(int64_t s, char **address, int *port, uint64_t desired_magic, int64_t justone);
#endif /* _RSOCKET_H_ */
