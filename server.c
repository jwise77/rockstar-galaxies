#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>
#include <time.h>
#include <inttypes.h>
#include <math.h>
#include <assert.h>
#include <sys/select.h>
#include <sys/wait.h>
#include "io/meta_io.h"
#include "io/io_bgc2.h"
#include "inet/rsocket.h"
#include "inet/address.h"
#include "config_vars.h"
#include "check_syscalls.h"
#include "fun_times.h"
#include "server.h"
#include "client.h"
#include "config.h"
#include "bounds.h"

float chunk_size[3];
int64_t chunks[3];
struct client_info *clients = NULL;
int64_t num_clients = 0;
int64_t time_start = 0;
int64_t server_error_state = 0;

#define for_writers(x) for (x=NUM_READERS; x<num_clients; x++)

void print_time(void) {
  int64_t time_now = time(NULL);
  fprintf(stderr, "[%6"PRId64"s] ", time_now-time_start);
}

void broadcast_msg(void *data, int64_t length) {
  for (int64_t i=num_clients-1; i>=0; i--)
    send_to_socket_delayconfirm(clients[i].cs, data, length);
  for (int64_t i=num_clients-1; i>=0; i--)
    send_to_socket_confirm(clients[i].cs, data, length);
}

void shutdown_clients() {
  for (int64_t i=0; i<num_clients; i++) {
    send_to_socket_noconfirm(clients[i].cs, "quit", 4);
    close_rsocket(clients[i].cs);
  }
}

void command_writers(char *cmd) {
  int64_t i;
  for_writers(i) send_to_socket_noconfirm(clients[i].cs, cmd, 4);
}

void send_to_writers(void *data, int64_t length) {
  int64_t i;
  for_writers(i) send_to_socket_noconfirm(clients[i].cs, data, length);
}

void accept_clients(int64_t s) {
  int port;
  int64_t c, num_writers=0, num_readers=0, accepted_client;
  uint64_t magic = ROCKSTAR_MAGIC;
  char cmd[5] = {0};

  timed_output("Accepting connections...\n");
  while (num_clients < NUM_READERS+NUM_WRITERS) {
    char *address=NULL;
    c = accept_connection(s, &address, &port);
    if (c < 0) continue;
    send_to_socket_noconfirm(c, &magic, sizeof(uint64_t));
    recv_from_socket(c, &magic, sizeof(uint64_t));

    if (magic != ROCKSTAR_MAGIC) {
      fprintf(stderr, "[Error] Received invalid client responses.  Check network connectivity.\n");
      shutdown_clients();
      exit(1);
    }

    recv_from_socket(c, cmd, 4);
    if ((num_readers < NUM_READERS) && 
	(!strcmp(cmd, "read") || !strcmp(cmd, "rdwr"))) {
      accepted_client = num_readers;
      num_readers++;
      send_to_socket_noconfirm(c, "read", 4);
    }
    else if ((num_writers < NUM_WRITERS) && 
	     (!strcmp(cmd, "writ") || !strcmp(cmd, "rdwr"))) {
      accepted_client = NUM_READERS + num_writers;
      send_to_socket_noconfirm(c, "writ", 4);
      num_writers++;
    }
    else {
      free(address);
      send_to_socket_noconfirm(c, "quit", 4);
      close_rsocket(c);
      continue;
    }
    clients[accepted_client].address = address;
    clients[accepted_client].port = port;
    clients[accepted_client].type = READER_TYPE;
    clients[accepted_client].cs = c;
    num_clients++;
  }
  timed_output("Accepted all reader / writer connections.\n");
}

int protocol_check(char *resp, char *expt) {
  if (!strcmp(resp, "err!")) {
    if (!server_error_state) {
      broadcast_msg("err!", 4);
      fprintf(stderr, "[Warning] Potentially fatal network error!\n");
    }
    server_error_state = 1;
    return 0;
  }
  else if (!strcmp(resp, "fail")) {
    fprintf(stderr, "[Error] Aborting because analysis process failed.\n");
    shutdown_clients();
    exit(1);
  }
  else if (strcmp(resp, expt)) {
    fprintf(stderr, "[Error] Protocol: expected %s, got %s!\n", expt, resp);
    shutdown_clients();
    exit(1);
  }
  return 1;
}

void wait_for_all_ready(int64_t client_min, int64_t client_max) {
  int64_t i;
  char cmd[5] = {0};
  for (i=client_min; i<client_max; i++)
    send_to_socket_delayconfirm(clients[i].cs, "rdy?", 4);
  for (i=client_min; i<client_max; i++) {
    send_to_socket_confirm(clients[i].cs, "rdy?", 4);
    recv_from_socket(clients[i].cs, cmd, 4);
    protocol_check(cmd, "rdy!");
    if (server_error_state) return;
  }
}

void reset_error() {
  int64_t i;
  char cmd[5] = {0};
  broadcast_msg("err!", 4);
  broadcast_msg("rdy?", 4);
  for (i=0; i<num_clients; i++) {
    if (recv_from_socket(clients[i].cs, cmd, 4) <= 0) {
      fprintf(stderr, "[Error] Fatal error: could not recover from network failure!\n");
      shutdown_clients();
      exit(1);
    }
    if (strcmp(cmd, "rdy!") != 0) i--;
  }
  server_error_state = 0;
}

void wait_for_done_signal() {
  char cmd[5] = {0};
  int64_t i;
  for_writers(i) {
    recv_from_socket(clients[i].cs, cmd, 4);
    protocol_check(cmd, "done");
  }
}

void command_writers_and_confirm(char *cmd) {
  command_writers(cmd);
  wait_for_done_signal();
}



void init_clients() {
  for (int64_t i=0; i<num_clients; i++)
    send_msg(clients[i].cs, clients[i].address, strlen(clients[i].address)+1);

  for (int64_t i=0; i<num_clients; i++)
    clients[i].serv_port = recv_msg_nolength(clients[i].cs, clients[i].serv_port);
  timed_output("Verified all reader / writer connections.\n");
}

int sort_by_address(const void *a, const void *b) {
  const struct client_info *c = a;
  const struct client_info *d = b;
  int res = strcmp(c->address, d->address);
  if (res) return res;
  return strcmp(c->serv_port, d->serv_port);
}

void sort_clients() {
  qsort(clients, NUM_READERS, sizeof(struct client_info), sort_by_address);
  qsort(clients+NUM_READERS, NUM_WRITERS, 
	sizeof(struct client_info), sort_by_address);
}

void transmit_client_info() {
  int64_t address_length = 0, i;
  for_writers(i) {
    address_length += strlen(clients[i].address)+1;
    address_length += strlen(clients[i].serv_port)+1;
  }
  char *client_addresses = check_realloc(NULL, sizeof(char)*address_length,
					 "Client addresses and ports");
  char *address_p = client_addresses;
  for_writers(i) {
    strcpy(address_p, clients[i].address);
    address_p += strlen(address_p)+1;
    strcpy(address_p, clients[i].serv_port);
    address_p += strlen(address_p)+1;
  }
  for_writers(i) {
    send_to_socket_noconfirm(clients[i].cs, "info", 4);
    send_to_socket(clients[i].cs, client_addresses, address_length);
  }
  free(client_addresses);
  timed_output("Transmitted all client connection info.\n");
}


void read_blocks(int64_t snap, int64_t pass) {
  int64_t block, reader, blocks_per_reader = NUM_BLOCKS / NUM_READERS;
  int64_t blocks_to_read = NUM_BLOCKS - NUM_READERS*pass; 
  if (NUM_BLOCKS % NUM_READERS) blocks_per_reader++;
  if (blocks_to_read > NUM_READERS) blocks_to_read = NUM_READERS;
  for (reader=0; reader < NUM_READERS; reader++) {
    block = reader*blocks_per_reader + pass;
    if (block >= NUM_BLOCKS) break;
    send_to_socket_noconfirm(clients[reader].cs, "snap", 4);
    send_to_socket_noconfirm(clients[reader].cs, &snap, sizeof(int64_t));
    send_to_socket_noconfirm(clients[reader].cs, "rdbk", 4);
    send_to_socket_noconfirm(clients[reader].cs, &block, sizeof(int64_t));
  }
  timed_output("Reading %"PRId64" blocks for snapshot %"PRId64"...\n",
	  blocks_to_read, snap);
}

#include "load_balance.c"

void decide_chunks() {
  factor_3(NUM_WRITERS, chunks);
  if (strlen(LOAD_BALANCE_SCRIPT)) decide_chunks_by_script();
  else if (NUM_WRITERS == 1) decide_chunks_for_volume_balance();
  else decide_chunks_for_memory_balance();
}

void decide_boundaries() {
  float box_size;
  int64_t reader;
  char cmd[5] = {0};
  for (reader=0; reader < NUM_READERS; reader++)
    send_to_socket_noconfirm(clients[reader].cs, "cnf?", 4);

  for (reader=0; reader < NUM_READERS; reader++) {
    recv_from_socket(clients[reader].cs, cmd, 4);
    protocol_check(cmd, "bxsz");
    recv_from_socket(clients[reader].cs, &box_size, sizeof(float));
    BOX_SIZE = box_size;
    recv_from_socket(clients[reader].cs, clients[reader].bounds,
		     sizeof(float)*6);
    recv_from_socket(clients[reader].cs, cmd, 4);
    protocol_check(cmd, "cnfg");
    recv_config(clients[reader].cs);
  }
  if ((BOX_SIZE < OVERLAP_LENGTH * 5) && PERIODIC) {
    shutdown_clients();
    fprintf(stderr, "[Error] Box size too small (%f) relative to overlap length (%f)!\n", BOX_SIZE, OVERLAP_LENGTH);
    exit(1);
  }
  wait_for_all_ready(0, NUM_READERS);
  decide_chunks();
}

void check_num_writers(void) {
  int64_t factors[3] = {0};
  factor_3(NUM_WRITERS, factors);
  if ((factors[0] < 2) || (factors[1] < 2) || (factors[2] < 2)) {
    fprintf(stderr, "[Error] NUM_WRITERS should be the product of at least three factors larger than 1 for periodic boundary conditions to be enabled!\n");
    fprintf(stderr, "[Error] (Currently, NUM_WRITERS = %"PRId64" = %"PRId64" x %"PRId64" x %"PRId64")\n", NUM_WRITERS, factors[0], factors[1], factors[2]);
    fprintf(stderr, "[Error] Please adjust NUM_WRITERS or set PERIODIC=0 in the config file.\n");
    exit(1);
  }
}

void send_bounds(int64_t i, int64_t j, float *bounds) {
  int64_t chunk = j-NUM_READERS;
  send_to_socket_noconfirm(clients[i].cs, bounds, sizeof(float)*6);
  send_to_socket_noconfirm(clients[i].cs, &chunk, sizeof(int64_t));
  send_to_socket_noconfirm(clients[i].cs, clients[j].address, strlen(clients[j].address)+1);
  send_to_socket_noconfirm(clients[i].cs, clients[j].serv_port, strlen(clients[j].serv_port)+1);
}

void transfer_data(int type, double overlap) {
  int64_t j;
  float bounds[6], *rcptbounds, *sndbounds;
  int64_t i = (type == DATA_PARTICLES) ? 0 : NUM_READERS;
  int64_t max_i = (type == DATA_PARTICLES) ? NUM_READERS : num_clients;
  char *send_cmd = (type == DATA_HALOS) ? "rcph" : "rcpt";
  if (type == DATA_BPARTICLES) send_cmd = "rcpb";
  if (type == DATA_HPARTICLES) send_cmd = "rchp";
  for (; i<max_i; i++) {
    for_writers(j) {
      rcptbounds = (type == DATA_HALOS) ? clients[j].bounds_prevsnap
	: clients[j].bounds;
      sndbounds = (type == DATA_HPARTICLES) ? clients[i].halo_bounds :
	clients[i].bounds;
      if (j==i && type == DATA_BPARTICLES) continue;
      if (bounds_overlap(sndbounds, rcptbounds, bounds, overlap)) {
	send_to_socket_noconfirm(clients[i].cs, send_cmd, 4);
	send_bounds(i,j,bounds);
      }
    }
  }
}

void transfer_particles() {
  int64_t i;
  char cmd[5] = {0};
  transfer_data(DATA_PARTICLES, 0);
  broadcast_msg("xfrp", 4);
  timed_output("Transferring particles to writers...\n");
  for (i=0; i<NUM_READERS; i++) {
    recv_from_socket(clients[i].cs, cmd, 4);
    protocol_check(cmd, "done");
    if (server_error_state) return;
  }
}

void transfer_bparticles() {
  transfer_data(DATA_BPARTICLES, 1.01*AVG_PARTICLE_SPACING*FOF_LINKING_LENGTH);
  timed_output("Transferring boundary particles between writers...\n");
  command_writers_and_confirm("xfrb");
}


void find_halos(int64_t snap) {
  int64_t i, chunk;
  for_writers(i) {
    send_to_socket_noconfirm(clients[i].cs, "done", 4);
    send_to_socket_noconfirm(clients[i].cs, "snap", 4);
    send_to_socket_noconfirm(clients[i].cs, &snap, sizeof(int64_t));
    send_to_socket_noconfirm(clients[i].cs, "cnfg", 4);
    send_config(clients[i].cs);
  }

  wait_for_all_ready(NUM_READERS, num_clients);
  if (server_error_state) return;

  if (!DUMP_PARTICLES[0]) {
    timed_output("Analyzing for FoF groups...\n");
    for_writers(i) {
      send_to_socket_noconfirm(clients[i].cs, "fofs", 4);
      chunk = i-NUM_READERS;
      send_to_socket_noconfirm(clients[i].cs, &chunk, sizeof(int64_t));
    }
    wait_for_done_signal();
    wait_for_all_ready(NUM_READERS, num_clients);
    transfer_bparticles();
    command_writers("done");
    wait_for_all_ready(NUM_READERS, num_clients);

    timed_output("Linking boundary particles...\n");
    command_writers_and_confirm("lnkb");
    command_writers("done");
    wait_for_all_ready(NUM_READERS, num_clients);

    timed_output("Analyzing for halos / subhalos...\n");
    command_writers("rock");
    load_balance();
    if (server_error_state) return;

    if (check_bgc2_snap(snap) || STRICT_SO_MASSES) {
      timed_output("Generating BGC2 files/SO Masses...\n");
      command_writers("hbnds");
      for_writers(i) recv_from_socket(clients[i].cs, clients[i].halo_bounds,
				      sizeof(float)*6);
      transfer_data(DATA_HPARTICLES, 0);
      command_writers_and_confirm("bgc2");
      command_writers("done");
    }
    command_writers("free");
  } else {
    timed_output("Dumping particles...\n");
    for_writers(i) {
      send_to_socket_noconfirm(clients[i].cs, "chnk", 4);
      chunk = i-NUM_READERS;
      send_to_socket_noconfirm(clients[i].cs, &chunk, sizeof(int64_t));
      send_to_socket_noconfirm(clients[i].cs, "outp", 4);
    }
    wait_for_done_signal();
  }
}

void get_bounds(int64_t snap)
{
  int64_t chunk = 0, i;
  for (chunk = 0; chunk<NUM_WRITERS; chunk++) {
    i = chunk+NUM_READERS;
    send_to_socket_noconfirm(clients[i].cs, "snap", 4);
    send_to_socket_noconfirm(clients[i].cs, &snap, sizeof(int64_t));
    send_to_socket_noconfirm(clients[i].cs, "chnk", 4);
    send_to_socket_noconfirm(clients[i].cs, &chunk, sizeof(int64_t));
    send_to_socket_noconfirm(clients[i].cs, "gbds", 4);
  }

  for (i = NUM_READERS; i<num_clients; i++)
    recv_from_socket(clients[i].cs, clients[i].bounds, sizeof(float)*6);
}

void _do_merger_tree_part2(int64_t snap) {
  int64_t timestep = 1, i;
  int64_t location = 0;
  char cmd[5] = {0};

  get_bounds(snap);
  for_writers(i) {
    send_to_socket_noconfirm(clients[i].cs, "rcph", 4);
    send_bounds(i,i,clients[i].bounds);
    send_to_socket_noconfirm(clients[i].cs, "xfrh", 4);
    send_to_socket_noconfirm(clients[i].cs, &timestep, sizeof(int64_t));
  }
  wait_for_done_signal();
  timed_output("Constructing merger tree...\n");
  for_writers(i) {
    send_to_socket_noconfirm(clients[i].cs, "done", 4);
    send_to_socket_noconfirm(clients[i].cs, "merg", 4);
  }
  command_writers("genc");
  for_writers(i) {
    recv_from_socket(clients[i].cs, &(clients[i].head_length), sizeof(int64_t));
    recv_from_socket(clients[i].cs, &(clients[i].cat_length), sizeof(int64_t));
    location += clients[i].head_length;
  }
  for_writers(i) {
    send_to_socket_noconfirm(clients[i].cs, "outc", 4);
    send_to_socket_noconfirm(clients[i].cs, &location, sizeof(int64_t));
    location += clients[i].cat_length;
    if (!PARALLEL_IO_CATALOGS) {
      send_to_socket_noconfirm(clients[i].cs, "rdy?", 4);
      recv_from_socket(clients[i].cs, cmd, 4);
      protocol_check(cmd, "rdy!");
    }
    send_to_socket_noconfirm(clients[i].cs, "delb", 4);
  }
  wait_for_all_ready(NUM_READERS, num_clients);
}

void do_merger_tree(int64_t snap) {
  int64_t i, timestep, starting_snap = 0;
  if ((SINGLE_SNAP && !snap) || (!SINGLE_SNAP && snap == STARTING_SNAP))
    starting_snap = 1;
  if (starting_snap && snap != NUM_SNAPS-1) return;

  timed_output("Loading merger tree information...\n");
  if (!starting_snap) { // Load in the current snapshot, with overlap
    timestep = 2;
    get_bounds(snap-1);
    for_writers(i)
      memcpy(clients[i].bounds_prevsnap, clients[i].bounds, sizeof(float)*6);

    get_bounds(snap);
    if (!p_bounds)
      check_realloc_s(p_bounds, sizeof(struct prev_bounds), NUM_WRITERS);
    for_writers(i)
      memcpy(p_bounds[i-NUM_READERS].bounds, clients[i].bounds,sizeof(float)*6);
    for_writers(i) {
      send_to_socket_noconfirm(clients[i].cs, "pbds", 4);
      send_to_socket_noconfirm(clients[i].cs, p_bounds,
			       sizeof(struct prev_bounds)*NUM_WRITERS);
    }

    transfer_data(DATA_HALOS, OVERLAP_LENGTH);
    if (server_error_state) return;

    for (i=NUM_READERS; i<num_clients; i++) {
      send_to_socket_noconfirm(clients[i].cs, "xfrh", 4);
      send_to_socket_noconfirm(clients[i].cs, &timestep, sizeof(int64_t));
    }
    wait_for_done_signal();
    command_writers("done");
    wait_for_all_ready(NUM_READERS, num_clients);

    if (server_error_state) return;
    _do_merger_tree_part2(snap-1);
  }
  if (snap == NUM_SNAPS-1) _do_merger_tree_part2(snap);
}


int64_t setup_server_port(void) {
  int64_t s, tries, addr_found = 0, auto_addr = 0, auto_port = 0;
#define s_address PARALLEL_IO_SERVER_ADDRESS
#define s_port PARALLEL_IO_SERVER_PORT
#define s_iface PARALLEL_IO_SERVER_INTERFACE
  if (!strcasecmp(s_address, "auto")) {
    auto_addr = 1;
    if (strlen(s_iface)) {
      s_address = get_interface_address(s_iface);
      if (s_address) addr_found = 1;
    }
    
    if (!addr_found) {
      s_address = check_realloc(NULL, 1024, "Allocating hostname.");
      if (gethostname(s_address, 1023)<0) {
	printf("Unable to get host address!\n");
	exit(1);
      }
    }
  }
  
  if (!strcasecmp(s_port, "auto")) {
    auto_port = 1;
    s_port = check_realloc(NULL, sizeof(char)*10, "Allocating port.");
    for (tries = 0; tries<500; tries++) {
      snprintf(s_port, 10, "%d", (rand()%63000)+2000);
      s = listen_at_addr(s_address, s_port);
      if (s>=0) break;
    }
  }
  else 
    s = listen_at_addr(s_address, s_port);

  if (s<0) {
    if (auto_addr || auto_port) { //Only server
      printf("Unable to start server on %s!\n", s_address);
      exit(1);
    }
    else return -1; //Must be a client
  }
  
  output_config("auto-rockstar.cfg");
  if (auto_addr) {
    free(s_address);
    s_address = "auto";
  }
  if (auto_port) {
    free(s_port);
    s_port = "auto";
  }

  return s;
}  

int server(void) {
  char buffer[1024];
  int64_t s, snap, num_passes, i, reload_parts = 0, n;
  int64_t data_size = sizeof(struct client_info)*(NUM_READERS+NUM_WRITERS);

  s = setup_server_port();
  if (s<0) return 0; //Client

  clients = check_realloc(clients,data_size, "Allocating client info.");
  memset(clients, 0, data_size);

  time_start = time(NULL);
  accept_clients(s);
  init_clients();
  sort_clients();
  transmit_client_info();
  num_passes = NUM_BLOCKS / NUM_READERS;
  if (NUM_BLOCKS % NUM_READERS) num_passes++;
  if (STARTING_SNAP > RESTART_SNAP) RESTART_SNAP = STARTING_SNAP;
  reload_parts = 1;

  for (snap = RESTART_SNAP; snap < NUM_SNAPS; snap++) {
    RESTART_SNAP = snap;
    output_config("restart.cfg");
    wait_for_all_ready(NUM_READERS, num_clients);
    if (!DO_MERGER_TREE_ONLY) {
      if (!PRELOAD_PARTICLES || reload_parts) {
	for (i=0; i<num_passes; i++) read_blocks(snap, i);
	reload_parts = 0;
      }
      decide_boundaries();
      transfer_particles();
      if (server_error_state) { reset_error(); reload_parts = 1; continue; }
      if (PRELOAD_PARTICLES && (snap < NUM_SNAPS-1)) 
	for (i=0; i<num_passes; i++) read_blocks(snap+1, i);
      find_halos(snap);
      if (server_error_state) { reset_error(); reload_parts = 1; continue; }
    }
    if (((strcasecmp(OUTPUT_FORMAT, "ASCII") != 0) || TEMPORAL_HALO_FINDING)
	&& !DUMP_PARTICLES[0] && !IGNORE_PARTICLE_IDS)
      do_merger_tree(snap);
    if (server_error_state) { reset_error(); reload_parts = 1; continue; }

    timed_output("[Success] Done with snapshot %"PRId64".\n", snap);
    if (strlen(RUN_PARALLEL_ON_SUCCESS)) {
      timed_output("Running external parallel analysis process for snapshot %"PRId64"...\n", snap);
      command_writers_and_confirm("rpos");
    }

    if (strlen(RUN_ON_SUCCESS)) {
      if (snapnames && snapnames[snap])
	snprintf(buffer, 1024, "%s %"PRId64" %s", 
		 RUN_ON_SUCCESS, snap, snapnames[snap]);
      else
	snprintf(buffer, 1024, "%s %"PRId64" %"PRId64, 
		 RUN_ON_SUCCESS, snap, snap);
      n = fork();
      if (n<=0) {
	if (system(buffer)!=0)
	  fprintf(stderr, "[Warning] Post-analysis command \"%s\" exited abnormally.\n", buffer);
      }
      if (n==0) exit(0);
    }
    if (SINGLE_SNAP) break;
  }
  while (wait(NULL)>=0);
  shutdown_clients();
  timed_output("[Finished]\n");
  return 1;
}
