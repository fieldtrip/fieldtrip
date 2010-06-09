#include <stdio.h>
#include <stdlib.h>
#include "peer.h"

/* this function returns 1 if the IP address corresponds with the local host
   or false if the IP address cannot be associated with the local host */

int localhost(const char *ipaddr)
{
		int family, s, status = 0, verbose = 0;
		struct ifaddrs *ifaddr, *ifa;
		char host[NI_MAXHOST];

		if (getifaddrs(&ifaddr) == -1) {
				perror("getifaddrs");
		}
		else {
				/* Walk through linked list, maintaining head pointer so we can free list later */

				for (ifa = ifaddr; ifa != NULL; ifa = ifa->ifa_next) {
						family = ifa->ifa_addr->sa_family;

						if (family == AF_INET || family == AF_INET6) {
								s = getnameinfo(ifa->ifa_addr,
												(family == AF_INET) ? sizeof(struct sockaddr_in) :
												sizeof(struct sockaddr_in6),
												host, NI_MAXHOST, NULL, 0, NI_NUMERICHOST);

								if (s == 0) {
										/* compare this hosts address with the user-specified address */
										status = (strcmp(host, ipaddr)==0);
										if (status)
												break; /* no reason to search further */
								}
						}
				} /* for looping over list */

				freeifaddrs(ifaddr);
		} /* if getifaddrs */

		if (verbose>0) {
				if (status)
						printf("<%s> == localhost\n", ipaddr);
				else
						printf("<%s> != localhost\n", ipaddr);
		}

		return status;
}

