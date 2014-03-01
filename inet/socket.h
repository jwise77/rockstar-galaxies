#ifndef _INET_SOCKET_H_
#define _INET_SOCKET_H_
#include <stdint.h>

#define RETRIES 10
#define TIMEOUT 2

void random_sleep(int timeout_secs);
int _connect_to_addr(char *host, char *port);
int _listen_at_addr(char *host, char *port);
int _accept_connection(int s, char **address, int *port);
int64_t _send_to_socket(int s, void *data, int64_t length);
int64_t _recv_from_socket(int s, void *data, int64_t length);
void *_recv_and_alloc(int s, void *data, int64_t length);
int64_t _send_msg(int s, void *data, int64_t length);
void *_recv_msg(int s, void *data, int64_t *length, int64_t offset);
void *_recv_msg_nolength(int s, void *data);
void set_network_io_error_cb(void (*cb)(int), int data);
void *socket_check_realloc(void *ptr, size_t size, char *reason);

#endif /* _INET_SOCKET_H_ */
