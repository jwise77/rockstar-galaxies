#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <string.h>
#include <sys/socket.h>
#include <sys/select.h>
#include <sys/poll.h>
#include <sys/signal.h>
#include <unistd.h>
#include <assert.h>
#include <errno.h>
#include "socket.h"
#include "rsocket.h"

//#define DEBUG_RSOCKET

#define SERVER_FLAG 1
#define NEW_FLAG 2
#define UNUSED_FLAG 4
#define RECEIVER_FLAG 8
#define SELECTED_FLAG 16
#define META_SELECTED_FLAG 32
#define DELAY_CONFIRM_FLAG 64

#define RPACKET_NO_CONFIRM_FLAG 1
#define RPACKET_DELAY_CONFIRM_FLAG 2
#define RPACKET_CONFIRM_SENT_FLAG 4

struct rsocket *rsockets = NULL;
int64_t num_rsockets = 0, poll_alloced = 0;
struct pollfd *rsocket_poll = NULL;

void rsocket_fatal(char *reason) {
  fprintf(stderr, "[Error] [Network] (PID: %d) %s\n", getpid(), reason);
  exit(1);
}

struct rsocket *rsocket_verify_id(int64_t id) {
  if (id > num_rsockets || id < 0)
    rsocket_fatal("Attempted use of invalid socket.");
  if (rsockets[id].flags & UNUSED_FLAG) {
    fprintf(stderr, "Socket id: %"PRId64"\n", id);
    fprintf(stderr, "Socket flags: %"PRId64"\n", rsockets[id].flags);
    rsocket_fatal("Attempted use of closed socket.");
  }
  return(rsockets+id);
}

void tag_rsocket(int64_t s) {
  struct rsocket *rs = rsocket_verify_id(s);
  rs->flags |= SELECTED_FLAG;
}

void clear_rsocket_tags(void) {
  for (int64_t i=0; i<num_rsockets; i++) {
    rsockets[i].flags -= (rsockets[i].flags & SELECTED_FLAG);
    rsockets[i].flags -= (rsockets[i].flags & META_SELECTED_FLAG);
  }
}

int64_t check_rsocket_tag(int64_t s) {
  if (s > num_rsockets || s < 0)
    rsocket_fatal("Attempted use of invalid socket.");
  if (rsockets[s].flags & UNUSED_FLAG) return 0;
  return (rsockets[s].flags & SELECTED_FLAG);
}

int64_t select_rsocket(int poll_type, double timeout) {
   int i, j, poll_events = 0, max_fd=0, max_i=0, nfds=0, res;

   if (poll_alloced != num_rsockets) {
     poll_alloced = num_rsockets;
     rsocket_poll = socket_check_realloc(rsocket_poll, sizeof(struct pollfd)*(num_rsockets+1),
					 "Allocating poll fds.");
   }
   if (poll_type & RSOCKET_READ) poll_events |= POLLIN;
   if (poll_type & RSOCKET_WRITE) poll_events |= POLLOUT;
   if (poll_type & RSOCKET_ERROR) poll_events |= POLLERR;
 
   //Check for new accepted sockets first
   if (poll_type & RSOCKET_READ) {
     for (i=0; i<num_rsockets; i++) {
       if ((rsockets[i].flags & NEW_FLAG) && !(rsockets[i].flags & UNUSED_FLAG)
	   && !(rsockets[rsockets[i].server_id].flags & UNUSED_FLAG) &&
	   (rsockets[rsockets[i].server_id].flags & SELECTED_FLAG) &&
	   (rsockets[rsockets[i].server_id].flags & SERVER_FLAG)) {
	 clear_rsocket_tags();
	 rsockets[rsockets[i].server_id].flags |= SELECTED_FLAG;
	 return(rsockets[i].server_id+1);
       }
     }
   }   
 
   //Flag server sockets to check for reaccepted connections:
   if (poll_type & RSOCKET_READ)
     for (i=0; i<num_rsockets; i++)
       if ((rsockets[i].flags & SELECTED_FLAG) && 
	   !(rsockets[i].flags & (UNUSED_FLAG)) &&
	   (rsockets[i].flags & RECEIVER_FLAG) && 
	   !(rsockets[rsockets[i].server_id].flags & UNUSED_FLAG))
	 rsockets[rsockets[i].server_id].flags |= META_SELECTED_FLAG;

   for (i=0; i<num_rsockets; i++) {
     if ((rsockets[i].flags & (SELECTED_FLAG|META_SELECTED_FLAG)) && 
	 !(rsockets[i].flags & (UNUSED_FLAG))) {
       if (rsockets[i].fd > max_fd) max_fd = rsockets[i].fd;
       rsocket_poll[nfds].fd = rsockets[i].fd;
       rsocket_poll[nfds].events = poll_events;
       nfds++;
     }
   }

   while ((res = poll(rsocket_poll, nfds, (timeout > 0.0) ? timeout*1000 : -1)) < 0
	  && (errno == EINTR));
   if (res < 0) {
     perror("[Warning] Socket poll() failed");
     clear_rsocket_tags();
     return 0;
   }
 
   for (j=0; j<nfds; j++)
     if (rsocket_poll[j].revents & (POLLERR|POLLNVAL))
       fprintf(stderr, "[Warning] error event mask %0x on fd %d\n", 
	       rsocket_poll[j].revents, rsocket_poll[j].fd);

   for (j=0; j<nfds; j++) {
     i = rsocket_from_fd(rsocket_poll[j].fd);
     assert(i>=0 && i<num_rsockets);
     if (rsocket_poll[j].revents & poll_events) {
       if ((rsockets[i].flags & SERVER_FLAG) &&
	   !(rsockets[i].flags & UNUSED_FLAG)) {
	 rsocket_accept_connection(i, NULL, NULL, 0, 1);
	 clear_rsocket_tags();
	 return 0;
       } else {
	 if (rsocket_poll[j].revents & POLLHUP) //Closed:
	   if (recv(rsockets[i].fd, &res, 1, MSG_PEEK)<=0) //Test for data
	     continue; //Skip if no data ready to read
	 rsockets[i].flags |= SELECTED_FLAG;
	 if (max_i < i) max_i = i;
       }
     } else {
       rsockets[i].flags -= (rsockets[i].flags & SELECTED_FLAG);
     }
   }
   return(max_i+1);
}

uint64_t gen_magic(void) {
  uint64_t magic = 0;
  while (!magic) {
    magic = ((uint64_t)rand()) << (uint64_t)32;
    magic += rand();
    for (int64_t i=0; i<num_rsockets; i++)
      if (rsockets[i].magic == magic) magic = 0;
  }
  return magic;
}

void close_rsocket(int64_t s) {
  if (s > num_rsockets || s < 0) 
    rsocket_fatal("Attempted closing of invalid socket.");
  if (rsockets[s].flags & UNUSED_FLAG) return;
  free(rsockets[s].address);
  free(rsockets[s].port);
  close(rsockets[s].fd);
  rsockets[s].flags = UNUSED_FLAG;
}

int rsocket_fd(int64_t s) {
  struct rsocket *rs = rsocket_verify_id(s);
  return rs->fd;
}

int64_t rsocket_from_fd(int fd) {
  for (int64_t i=0; i<num_rsockets; i++)
    if (rsockets[i].fd == fd && !(rsockets[i].flags & UNUSED_FLAG))
      return i;
  return -1;
}

int64_t add_rsocket(struct rsocket s) {
  int64_t u = -1;
  for (s.id=0; s.id<num_rsockets; s.id++) {
    if (rsockets[s.id].flags & UNUSED_FLAG) u=s.id;
    if (s.magic &&
	(rsockets[s.id].magic == s.magic) &&
	(rsockets[s.id].flags & RECEIVER_FLAG)) {
      if (!(rsockets[s.id].flags & UNUSED_FLAG)) {
	if (rsockets[s.id].address != s.address) free(rsockets[s.id].address);
	if (rsockets[s.id].port != s.port) free(rsockets[s.id].port);
	s.sseq=rsockets[s.id].sseq;
	s.rseq=rsockets[s.id].rseq;
      }
      rsockets[s.id] = s;
      return s.id;
    }
  }

  if (u>=0) {
    rsockets[u] = s;
    rsockets[u].id = u;
    return u;
  }
  rsockets = socket_check_realloc(rsockets, sizeof(struct rsocket)*(num_rsockets+1), "Adding new rsocket.");
  s.id = num_rsockets;
  rsockets[num_rsockets] = s;
  num_rsockets++;
  return (s.id);
}

void _reconnect_to_addr(struct rsocket *n) {
  uint64_t magic = RSOCKET_VERIFY;
  n->fd = _connect_to_addr(n->address, n->port); //Always succeeds

  if (!_send_to_socket(n->fd, &(magic), sizeof(uint64_t)) || 
      (!_send_to_socket(n->fd, &(n->magic), sizeof(uint64_t))))
    rsocket_fatal("Couldn't send magic number over socket.");
  if (!_recv_from_socket(n->fd, &magic, sizeof(uint64_t)))
    rsocket_fatal("Couldn't receive magic number over socket.");
  if (!n->magic) n->magic = magic;
  else if (n->magic != magic) rsocket_fatal("Couldn't reconnect to address!");  
}

int64_t connect_to_addr(char *host, char *port) {
  struct rsocket n = {0};
  if (!num_rsockets) signal(SIGPIPE, SIG_IGN);
  int64_t c = add_rsocket(n);
  rsockets[c].address = strdup(host);
  rsockets[c].port = strdup(port);
  _reconnect_to_addr(rsockets+c);
  return c;
}

int64_t listen_at_addr(char *host, char *port) {
  if (!num_rsockets) signal(SIGPIPE, SIG_IGN);
  struct rsocket ns = {0};
  ns.fd = _listen_at_addr(host, port);
  if (ns.fd<0) return ns.fd;
  ns.magic = gen_magic();
  ns.flags = SERVER_FLAG;
  ns.address = strdup(host);
  ns.port = strdup(port);
  int64_t ret = add_rsocket(ns);
  assert(ns.address == rsockets[ret].address &&
	 ns.port == rsockets[ret].port);
  return ret;
}

int64_t rsocket_accept_connection(int64_t s, char **address, int *port, uint64_t desired_magic, int64_t justone) {
  struct rsocket *ts = rsocket_verify_id(s);
  char port_buffer[20] = {0};
  int our_port;
  uint64_t magic;
  if (!(ts->flags & SERVER_FLAG))  rsocket_fatal("Not a server socket!");
  while (1) {
    struct rsocket nc = {0};
    ts = rsocket_verify_id(s);
    nc.fd = _accept_connection(ts->fd, &(nc.address), &our_port);
    if ((_recv_from_socket(nc.fd, &(magic), sizeof(uint64_t)) != sizeof(uint64_t))
	|| (magic != RSOCKET_VERIFY) ||
	(_recv_from_socket(nc.fd, &(nc.magic), sizeof(uint64_t)) != sizeof(uint64_t))){
      fprintf(stderr, "[Warning] Ignoring non-rsocket client.\n");
      close(nc.fd); //Ignore connections from non-rsocket clients
      free(nc.address);
      if (justone) return -1;
      continue;
    }

    nc.flags |= RECEIVER_FLAG;
    if (address && (address != &(nc.address))) *address = strdup(nc.address);
    if (port) *port = our_port;
    snprintf(port_buffer, 20, "%d", our_port);
    nc.port = strdup(port_buffer);
    nc.server_id = s;
    int64_t id = add_rsocket(nc);
#ifdef DEBUG_RSOCKET
    fprintf(stderr, "Accepted connection %"PRId64" with sseq %"PRId64" rseq %"PRId64"\n",
	    id, rsockets[id].sseq, rsockets[id].rseq);
#endif /*DEBUG_RSOCKET */
    if (!rsockets[id].magic) rsockets[id].magic = gen_magic();
    _send_to_socket(nc.fd, &(rsockets[id].magic), sizeof(uint64_t));
    if (justone && !nc.magic) rsockets[id].flags |= NEW_FLAG;
    if (nc.magic==desired_magic || justone) return id;
    if (!nc.magic) rsockets[id].flags |= NEW_FLAG;
  }
}

int64_t accept_connection(int64_t s, char **address, int *port) {
  for (int64_t i=0; i<num_rsockets; i++) {
    if (rsockets[i].flags & NEW_FLAG) {
      rsockets[i].flags -= (rsockets[i].flags & NEW_FLAG);
      return i;
    }
  }
  return(rsocket_accept_connection(s, address, port, 0, 0));
}


struct rsocket *repair_connection(int64_t *s) {
  struct rsocket *rs = rsocket_verify_id(*s);
  close(rs->fd);
  if (rs->flags & RECEIVER_FLAG)
    *s = rsocket_accept_connection(rs->server_id, NULL, NULL, rs->magic, 0);
  else _reconnect_to_addr(rsockets + *s);
  return rsocket_verify_id(*s);
}


int64_t send_to_rsocket(int64_t s, void *data, int64_t length, int64_t flags) {
  int64_t confirm=0;
  struct rsocket *rs = rsocket_verify_id(s);
  struct rpacket_header rp;
  int64_t retry = 0;
  rp.length = length;
  rp.seq = rs->sseq;
  rp.flags = flags;

  if ((rs->flags & DELAY_CONFIRM_FLAG) && !(flags & RPACKET_CONFIRM_SENT_FLAG))
    rsocket_fatal("Programmer broke promise to confirm receipt of packet!");
  if ((flags & RPACKET_CONFIRM_SENT_FLAG) && !(rs->flags & DELAY_CONFIRM_FLAG))
    flags -= (flags & RPACKET_CONFIRM_SENT_FLAG);

  while (1) {
    retry++;

    if (((flags & RPACKET_CONFIRM_SENT_FLAG) ||
	 ((_send_to_socket(rs->fd, &(rs->magic), sizeof(uint64_t))==sizeof(uint64_t))
	  && (_send_to_socket(rs->fd, &rp, sizeof(struct rpacket_header))==
	  sizeof(struct rpacket_header)) &&
	  (_send_to_socket(rs->fd, data, length)==length))) &&
	((flags & (RPACKET_NO_CONFIRM_FLAG | RPACKET_DELAY_CONFIRM_FLAG)) ||
	 ((_recv_from_socket(rs->fd, &confirm, sizeof(int64_t))==sizeof(int64_t)) &&
	  (confirm==1)))) {
      if (!(flags & RPACKET_DELAY_CONFIRM_FLAG)) {
#ifdef DEBUG_RSOCKET
	fprintf(stderr, "[Sent] %"PRId64" bytes to socket %"PRId64", seq %"PRId64", seqnow %"PRId64". ", rp.length, s, rp.seq, rs->sseq);
	if (length < 100) fwrite(data, 1, rp.length, stderr);
	fprintf(stderr, "\n");
#endif /*DEBUG_RSOCKET*/
	rs->sseq++;
	rs->flags -= (rs->flags & DELAY_CONFIRM_FLAG);
      }
      else rs->flags |= DELAY_CONFIRM_FLAG;
      break;
    }
    flags -= (flags & RPACKET_CONFIRM_SENT_FLAG);
    fprintf(stderr, "[Network] Packet send retry count at: %"PRId64"\n", retry);
    rs = repair_connection(&s); 
  }
  rs->last_data = data;
  return length;
}

int64_t send_to_socket(int64_t s, void *data, int64_t length) {
  return send_to_rsocket(s, data, length, 0);
}

int64_t send_to_socket_noconfirm(int64_t s, void *data, int64_t length) {
  return send_to_rsocket(s, data, length, RPACKET_NO_CONFIRM_FLAG);
}

int64_t send_to_socket_delayconfirm(int64_t s, void *data, int64_t length) {
  return send_to_rsocket(s, data, length, RPACKET_DELAY_CONFIRM_FLAG);
}

int64_t send_to_socket_confirm(int64_t s, void *data, int64_t length) {
  return send_to_rsocket(s, data, length, RPACKET_CONFIRM_SENT_FLAG);
}


int64_t recv_from_rsocket(int64_t s, void *data, int64_t length, int64_t offset, int64_t alloc) {
  int64_t recvd, lcount=0, to_read = 0, confirm=1;
  uint64_t magic;
  char buffer[1600];
  struct rsocket *rs = rsocket_verify_id(s);
  struct rpacket_header rp;
  int64_t retry = 0;

  if (!length && !alloc) return 0;
  if ((rs->flags & DELAY_CONFIRM_FLAG))
    rsocket_fatal("Programmer broke promise to confirm receipt of packet!");

  while (1) {
    retry++;
    if ((_recv_from_socket(rs->fd, &magic, sizeof(uint64_t))==sizeof(uint64_t)) &&
	(magic == rs->magic) &&
	(_recv_from_socket(rs->fd, &rp, sizeof(struct rpacket_header))
	 ==sizeof(struct rpacket_header))) {
      //Check for duplicate sequence:
      if (rp.seq < rs->rseq && rs->rseq < (INT64_MAX - 10)) {
	fprintf(stderr, "[Warning] Ignoring duplicate sequence %"PRId64" (seqnow: %"PRId64")\n",
		rp.seq, rs->rseq);
	while (lcount < rp.length) {
	  to_read = rp.length - lcount;
	  if (to_read > 1600) to_read = 1600;
	  recvd = _recv_from_socket(rs->fd, buffer, to_read);
	  if (recvd != to_read) {
	    fprintf(stderr, "[Warning] Incomplete packet received!\n");
	    break;
	  }
	  lcount += recvd;
	}
	if (lcount == rp.length) {
	  if (!(rp.flags & RPACKET_NO_CONFIRM_FLAG))
	    _send_to_socket(rs->fd, &confirm, sizeof(int64_t));
	  continue;
	}
      }
      else {
	if (alloc) {
	  if (rp.length+offset > length)  {
	    length = rp.length+offset;
	    data = socket_check_realloc(data, length, "receive buffer");
	  }
	}
	else if ((length != rp.length+offset)) {
	  fprintf(stderr, "Expected receive length (%"PRId64") != actual receive length (%"PRId64")\n", length, rp.length+offset);
	  if (length < 100) {
	    fprintf(stderr, "Actual data received: ");
	    for (int64_t i=0; i<length; i++) fprintf(stderr, "%c", ((char *)data)[i]);
	    fprintf(stderr, "\n");
	    assert(0);
	  }
	  exit(1);
	}
	if (_recv_from_socket(rs->fd, ((char *)data)+offset, rp.length)==rp.length) {
	   if (!(rp.flags & RPACKET_NO_CONFIRM_FLAG))
	     _send_to_socket(rs->fd, &confirm, sizeof(int64_t));
	   rs->last_data = data;
#ifdef DEBUG_RSOCKET
	   fprintf(stderr, "[Received] %"PRId64" bytes from socket %"PRId64", seq %"PRId64", seqnow %"PRId64". ", rp.length, s, rp.seq, rs->rseq);
	   if (length < 100) fwrite(((char *)data)+offset, 1, rp.length, stderr);
	   fprintf(stderr, "\n");
#endif /*DEBUG_RSOCKET*/
	   rs->rseq = rp.seq+1;
	   break;
	}
      }
    }
    fprintf(stderr, "[Network] Packet receive retry count at: %"PRId64"\n", retry);
    rs = repair_connection(&s);
  }
  return length;
}

int64_t send_msg(int64_t s, void *data, int64_t length) {
  return send_to_socket(s, data, length);
}

int64_t recv_from_socket(int64_t s, void *data, int64_t length) {
  return recv_from_rsocket(s, data, length, 0, 0);
}

void *recv_and_alloc(int64_t s, void *data, int64_t length) {
  data = socket_check_realloc(data, length, "receive buffer");
  recv_from_rsocket(s, data, length, 0, 0);
  return data;
}

void *recv_msg(int64_t s, void *data, int64_t *length, int64_t offset) {
  *length = recv_from_rsocket(s, data, *length, offset, 1);
  return rsockets[s].last_data;
}

void *recv_msg_nolength(int s, void *data) {
  recv_from_rsocket(s, data, 0, 0, 1);
  return rsockets[s].last_data;
}
