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

#include <sys/types.h>
#include <time.h>
#include <string.h>
#include <stdio.h>

#include "peer.h"
#include "extern.h"
#include "platform_includes.h"

#ifdef COMPILER_MINGW
/* TODO: replace inet_ntop by getnameinfo(sa, salen, ipaddr, maxlen, NULL, 0, 0)
   for all platforms
*/
char *inet_ntop(int af, const void *src, char *dst, socklen_t maxlen) {
	if (af != AF_INET || maxlen < INET_ADDRSTRLEN) {
		dst[0] = 0;
		return NULL;
	} else {
		unsigned char *addr = (unsigned char *) src;
		snprintf(dst, maxlen, "%u.%u.%u.%u", addr[0], addr[1], addr[2], addr[3]);
		return dst;
	}
}
#endif

typedef struct {
		void **discovery;
		int *fd;
} threadlocal_t;

void cleanup_discover(void *arg) {
		peerlist_t *next = NULL;
		threadlocal_t *threadlocal = NULL;
		threadlocal = (threadlocal_t *)arg;

		DEBUG(LOG_DEBUG, "cleanup_discover()");

		if (discoverStatus==0)
				return;

		if (threadlocal && *threadlocal->discovery) {
				FREE(*threadlocal->discovery);
		}

		if (threadlocal && (*threadlocal->fd)>0) {
				closesocket(*threadlocal->fd);
				*threadlocal->fd = 0;
		}

		pthread_mutex_lock(&mutexpeerlist);
		while (peerlist) {
				next = peerlist->next;
				FREE(peerlist->host);
				FREE(peerlist);
				peerlist = next;
		}
		pthread_mutex_unlock(&mutexpeerlist);

		pthread_mutex_lock(&mutexstatus);
		discoverStatus = 0;
		pthread_mutex_unlock(&mutexstatus);

		pthread_mutex_lock(&mutexthreadcount);
		threadcount--;
		pthread_mutex_unlock(&mutexthreadcount);

		return;
}

void *discover(void *arg) {
		int accept = 1;
		int fd = 0;
		int i = 0;
		int localhost = 0;
		int nbytes, found = 0;
		int one = 1;
		hostdef_t  *discovery = NULL;
		peerlist_t *peer = NULL, *next = NULL;
		unsigned int addrlen;

		/* these variables are for the socket */
		struct sockaddr_in addr;
		struct ip_mreq mreq;
		int optval;

		threadlocal_t threadlocal;
		threadlocal.discovery = &discovery;
		threadlocal.fd        = &fd;

		/* this is for debugging */
		pthread_mutex_lock(&mutexthreadcount);
		threadcount++;
		pthread_mutex_unlock(&mutexthreadcount);

		pthread_cleanup_push(cleanup_discover, &threadlocal);

		/* the status contains the thread id when running, or zero when not running */
		pthread_mutex_lock(&mutexstatus);
		if (discoverStatus==0) {
				discoverStatus = 1;
				/* signal that this thread has started */
				pthread_cond_signal(&condstatus);
				pthread_mutex_unlock(&mutexstatus);
		}
		else {
				pthread_mutex_unlock(&mutexstatus);
				goto cleanup;
		}

		/* create what looks like an ordinary UDP socket */
		if ((fd=socket(AF_INET,SOCK_DGRAM,0)) < 0) {
				perror("discover socket");
				DEBUG(LOG_ERR, "error: discover socket");
				goto cleanup;
		}

		/* allow multiple sockets to use the same PORT number */
		if (setsockopt(fd,SOL_SOCKET,SO_REUSEPORT,&one,sizeof(one)) < 0) {
				perror("discover setsockopt so_reuseaddr");
				DEBUG(LOG_ERR, "error: discover setsockopt so_reuseaddr");
				goto cleanup;
		}

		/* set up destination address */
		memset(&addr,0,sizeof(addr));
		addr.sin_family=AF_INET;
		addr.sin_addr.s_addr=htonl(INADDR_ANY); /* N.B.: differs from sender */
		addr.sin_port=htons(ANNOUNCE_PORT);

		/* prevent "bind: Address already in use" */
		optval = 1;
		if (setsockopt(fd, SOL_SOCKET, SO_REUSEADDR, (const char*)&optval, sizeof(optval)) < 0) {
				perror("discover setsockopt");
				DEBUG(LOG_ERR, "error: discover setsockopt");
				goto cleanup;
		}

		/* bind to receive address */
		if (bind(fd,(struct sockaddr *)&addr,sizeof(addr)) < 0) {
				perror("discover bind");
				DEBUG(LOG_ERR, "error: discover bind");
				goto cleanup;
		}

		/* use setsockopt() to request that the kernel join a multicast group */
		mreq.imr_multiaddr.s_addr=inet_addr(ANNOUNCE_GROUP);
		mreq.imr_interface.s_addr=htonl(INADDR_ANY);
		if (setsockopt(fd,IPPROTO_IP,IP_ADD_MEMBERSHIP,&mreq,sizeof(mreq)) < 0) {
				perror("discover setsockopt ip_add_membership");
				DEBUG(LOG_ERR, "error: discover setsockopt ip_add_membership");
				goto cleanup;
		}

		/* now just enter an infinite loop */
		while (1) {
				char ipaddr[INET_ADDRSTRLEN];
				
				if ((discovery = malloc(sizeof(hostdef_t)))==NULL) {
						perror("discover malloc");
						DEBUG(LOG_ERR, "error: discover malloc");
						goto cleanup;
				}

				/* note that this might be thread cancelation point, but I am not sure */
				addrlen=sizeof(addr);
				if ((nbytes=recvfrom(fd,discovery,sizeof(hostdef_t),0,(struct sockaddr *)&addr,&addrlen)) < 0) {
						perror("discover recvfrom");
						DEBUG(LOG_ERR, "error: discover recvfrom");
						goto cleanup;
				}

				if (discovery->version!=VERSION) {
						FREE(discovery);
						continue;
				}

				/* the UDP packet specifies the IP address of the sender */
				
				inet_ntop(AF_INET, &addr.sin_addr, ipaddr, INET_ADDRSTRLEN);

				/* there seems to be a thread cancelation point inside check_localhost  */
				/* therefore it should be called only when all mutexes are unlocked */
				localhost = check_localhost(ipaddr);

				/*
				   DEBUG(LOG_DEBUG, "discover: host->name = %s", discovery->name);
				   DEBUG(LOG_DEBUG, "discover: host->port = %u", discovery->port);
				   DEBUG(LOG_DEBUG, "discover: host->id   = %u", discovery->id);
				   DEBUG(LOG_DEBUG, "discover: IP address = %s", ipaddr);
				 */

				/* check whether the peer should be listed */
				accept = 1;
				accept = (accept && ismember_userlist (discovery->user));
				accept = (accept && ismember_grouplist(discovery->group));
				accept = (accept && ismember_hostlist (discovery->name));

				if (!accept) {
						DEBUG(LOG_DEBUG, "discover: rejecting %s:%u, id = lu", discovery->name, discovery->port, discovery->id);
						FREE(discovery);
						continue;
				}
				else {
						DEBUG(LOG_DEBUG, "discover: accepting %s:%u, id = %u", discovery->name, discovery->port, discovery->id);
				}

				pthread_mutex_lock(&mutexpeerlist);

				/* remove previous observations of this discovery from the list */
				/* test the first item on the list */
				if (peerlist)	{
						found = (peerlist->host->id==discovery->id);
						if (found) {
								/* delete the first item in the list */
								next = peerlist->next;
								FREE(peerlist->host);
								FREE(peerlist);
								peerlist = next;
						}
				}

				/* remove previous observations of this discovery from the list */
				/* traverse the list */
				peer = peerlist;
				while(peer) {
						/* test the next item on the list */
						next = peer->next;
						if (next) {
								found = (next->host->id==discovery->id);
								if (found) {
										/* delete the next item in the list */
										peer->next = next->next;
										FREE(next->host);
										FREE(next);
										break;
								}
						}
						peer = peer->next;
				}

				/* add the new discovery to the list */
				peer       = (peerlist_t *)malloc(sizeof(peerlist_t));
				peer->host = (hostdef_t *)malloc(sizeof(hostdef_t));
				memcpy(peer->host, discovery, sizeof(hostdef_t));

				FREE(discovery);

				if (localhost)
						/* use the loopback IP address instead of the external IP address */
						strncpy(peer->ipaddr, "127.0.0.1", INET_ADDRSTRLEN);
				else
						/* use the external IP address */
						strncpy(peer->ipaddr, ipaddr, INET_ADDRSTRLEN);

				peer->time      = time(NULL);
				peer->next      = peerlist;
				peerlist        = peer;

				/* give some debug output
				   i = 0;
				   peer = peerlist;
				   while(peer) {
				   DEBUG(LOG_DEBUG, "discover: peerlist[%d] =", i);
				   DEBUG(LOG_DEBUG, "discover:   host.name = %s", peer->host->name);
				   DEBUG(LOG_DEBUG, "discover:   host.port = %u", peer->host->port);
				   DEBUG(LOG_DEBUG, "discover:   host.id   = %u", peer->host->id);
				   DEBUG(LOG_DEBUG, "discover:   ipaddr    = %s", peer->ipaddr);
				   DEBUG(LOG_DEBUG, "discover:   time      = %s", ctime(&(peer->time)));
				   peer = peer->next ;
				   i++;
				   }
				 */

				pthread_mutex_unlock(&mutexpeerlist);

				/* note that this is a thread cancelation point */
				pthread_testcancel();

		} /* while (1) */

cleanup:
		printf(""); /* otherwise the pthread_cleanup_pop won't compile */
		pthread_cleanup_pop(1);
		return NULL;
}
