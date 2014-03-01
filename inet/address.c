#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <arpa/inet.h>
#include <sys/socket.h>
#include <netdb.h>
#include <ifaddrs.h>
#include <unistd.h>

#ifndef   NI_MAXHOST
#define   NI_MAXHOST 1025
#endif

//Adapted from the getifaddrs(7) man page
char *get_interface_address(char *ifname) {
  struct ifaddrs *ifaddr, *ifa;
  char *address = 0, host[NI_MAXHOST];
  int family, s;

  s = getifaddrs(&ifaddr);

  if (s != -1) {
    for (ifa = ifaddr; ifa != NULL; ifa = ifa->ifa_next)
      {
	if (ifa->ifa_addr == NULL) continue;
	if (strcmp(ifa->ifa_name, ifname)) continue;
	family = ifa->ifa_addr->sa_family;
	if (family == AF_INET || family == AF_INET6) {
	  s = getnameinfo(ifa->ifa_addr,
			  (family == AF_INET) ? sizeof(struct sockaddr_in) :
			  sizeof(struct sockaddr_in6),
			  host, NI_MAXHOST, NULL, 0, NI_NUMERICHOST);
	  if (s == 0) {
	    fprintf(stderr, "[Note   ] Using %s for server address on iface %s.\n", host, ifname);
	    address = strdup(host);
	    break;
	  }
	}
      }
    freeifaddrs(ifaddr);
  }

  if (!address)
    fprintf(stderr, "[Warning] Unable to get interface addresses; using hostname instead.\n");
  return address;
}
