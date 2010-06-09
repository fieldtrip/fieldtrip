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

typedef struct {
		void *discovery;
		int fd;
} threadlocal_t;

void cleanup_discover(void *arg) {
		threadlocal_t *threadlocal = NULL;
		peerlist_t *next = NULL;

		threadlocal = (threadlocal_t *)arg;
		if (threadlocal && threadlocal->discovery) {
				FREE(threadlocal->discovery);
		}
		if (threadlocal && threadlocal->fd>0) {
				close(threadlocal->fd);
				threadlocal->fd = -1;
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
		int i = 0;
		int fd = 0;
		unsigned int addrlen;
		int nbytes, verbose = 0, found = 0;
		int one = 1;
		int accept = 1;
		peerlist_t *peer = NULL, *next = NULL;
		hostdef_t  *discovery = NULL;

		/* these variables are for the socket */
		struct sockaddr_in addr;
		struct ip_mreq mreq;
		int optval;

		threadlocal_t threadlocal;
		threadlocal.discovery = NULL;
		threadlocal.fd = -1;

		/* this is for debugging */
		pthread_mutex_lock(&mutexthreadcount);
		threadcount++;
		pthread_mutex_unlock(&mutexthreadcount);

		pthread_cleanup_push(cleanup_discover, NULL);

		/* the status contains the thread id when running, or zero when not running */
		pthread_mutex_lock(&mutexstatus);
		if (discoverStatus==0) {
				discoverStatus = 1;
				pthread_mutex_unlock(&mutexstatus);
		}
		else {
				pthread_mutex_unlock(&mutexstatus);
				goto cleanup;
		}

		/* create what looks like an ordinary UDP socket */
		if ((fd=socket(AF_INET,SOCK_DGRAM,0)) < 0) {
				perror("discover socket");
				goto cleanup;
		}

		/* this will be closed at cleanup */
		threadlocal.fd = fd;

		/* allow multiple sockets to use the same PORT number */
		if (setsockopt(fd,SOL_SOCKET,SO_REUSEPORT,&one,sizeof(one)) < 0) {
				perror("discover setsockopt so_reuseaddr");
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
				goto cleanup;
		}

		/* bind to receive address */
		if (bind(fd,(struct sockaddr *)&addr,sizeof(addr)) < 0) {
				perror("discover bind");
				goto cleanup;
		}

		/* use setsockopt() to request that the kernel join a multicast group */
		mreq.imr_multiaddr.s_addr=inet_addr(ANNOUNCE_GROUP);
		mreq.imr_interface.s_addr=htonl(INADDR_ANY);
		if (setsockopt(fd,IPPROTO_IP,IP_ADD_MEMBERSHIP,&mreq,sizeof(mreq)) < 0) {
				perror("discover setsockopt ip_add_membership");
				goto cleanup;
		}

		/* now just enter an infinite loop */
		while (1) {

				addrlen=sizeof(addr);

				if ((discovery = malloc(sizeof(hostdef_t)))==NULL) {
						perror("discover malloc");
						goto cleanup;
				}

				/* this will be deallocated at cleanup */
				threadlocal.discovery = discovery;

				/* note that this might be thread cancelation point, but I am not sure */
				if ((nbytes=recvfrom(fd,discovery,sizeof(hostdef_t),0,(struct sockaddr *)&addr,&addrlen)) < 0) {
						perror("discover recvfrom");
						goto cleanup;
				}

				if (discovery->version!=VERSION) {
						FREE(discovery);
						continue;
				}

				if (verbose>0) {
						fprintf(stderr, "\n");
						fprintf(stderr, "discover: host->name = %s\n", discovery->name);
						fprintf(stderr, "discover: host->port = %d\n", discovery->port);
						fprintf(stderr, "discover: host->id   = %d\n", discovery->id);
				}

				/* check whether the peer should be listed */
				accept = 1;
				accept = (accept & ismember_userlist (discovery->user));
				accept = (accept & ismember_grouplist(discovery->group));
				accept = (accept & ismember_hostlist (discovery->name));

				if (!accept) {
						FREE(discovery);
						continue;
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

				/* the UDP packet specifies the IP address of the sender */
				inet_ntop(AF_INET, &addr.sin_addr, peer->ipaddr, INET_ADDRSTRLEN);
				/* if possible use the loopback IP address instead of external IP address */
				if (localhost(peer->ipaddr)==1)
					strncpy(peer->ipaddr, "127.0.0.1", INET_ADDRSTRLEN);
				peer->time      = time(NULL);
				peer->next      = peerlist;
				peerlist        = peer;

				if (verbose>1) {
						i = 0;
						peer = peerlist;
						while(peer) {
								fprintf(stderr, "discover: peerlist[%d] =\n", i);
								fprintf(stderr, "discover:   host.name = %s\n", peer->host->name);
								fprintf(stderr, "discover:   host.port = %d\n", peer->host->port);
								fprintf(stderr, "discover:   host.id   = %d\n", peer->host->id);
								fprintf(stderr, "discover:   time      = %s",   ctime(&(peer->time)));
								peer = peer->next ;
								i++;
						}
				}

				pthread_mutex_unlock(&mutexpeerlist);

				/* note that this is a thread cancelation point */
				pthread_testcancel();


		} /* while (1) */

cleanup:
		printf(""); /* otherwise the pthread_cleanup_pop won't compile */
		pthread_cleanup_pop(1);
		return NULL;
}
