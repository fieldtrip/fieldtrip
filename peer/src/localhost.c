/*
 * Copyright (C) 2010, Robert Oostenveld
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/
 *
 */

#include <stdio.h>
#include <stdlib.h>

#include "platform_includes.h"
#include "peer.h"
#include "extern.h"

/* this function returns 1 if the IP address corresponds with the local host
   or 0 if the IP address cannot be associated with the local host */
int check_localhost(const char *ipaddr)
{
#if defined (PLATFORM_WIN32) || defined(PLATFORM_WIN64)
		return 0;

#elif defined (PLATFORM_LINUX) || defined(PLATFORM_OSX)
		int family, s, found = 0;
		struct ifaddrs *ifaddr = NULL;
		struct ifaddrs *ifa;
		char host[NI_MAXHOST];

		/* get the list with network interfaces */
		if (getifaddrs(&ifaddr) == -1) {
				perror("getifaddrs");
				DEBUG(LOG_ERR, "error: getifaddrs");
		}

		/* walk through the linked list, maintaining head pointer so we can free list later */
		for (ifa = ifaddr; ifa != NULL; ifa = ifa->ifa_next) {
				/* the following line fails when -fpack-struct is used during compilation */
				family = ifa->ifa_addr->sa_family;

				if (family == AF_INET)
						s = getnameinfo(ifa->ifa_addr, sizeof(struct sockaddr_in ), host, NI_MAXHOST, NULL, 0, NI_NUMERICHOST);
				else if (family == AF_INET6)
						s = getnameinfo(ifa->ifa_addr, sizeof(struct sockaddr_in6), host, NI_MAXHOST, NULL, 0, NI_NUMERICHOST);
				else
						s = -1;

				/* compare this hosts address with the user-specified address */
				found = (s==0) && (strcmp(host, ipaddr)==0);

				if (found)
						break;
		} /* for looping over list */

		freeifaddrs(ifaddr);

		if (found)
				DEBUG(LOG_DEBUG, "localhost: <%s>", ipaddr);

		return found;
#endif
}

