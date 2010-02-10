#include <sys/types.h>
#include <time.h>
#include <string.h>
#include <stdio.h>

#include "peer.h"
#include "extern.h"
#include "platform_includes.h"

typedef struct {
		void *host;
		int fd;
} threadlocal_t;

void cleanup_discover(void *arg) {
		threadlocal_t *threadlocal = NULL;
		peerlist_t *next = NULL;

    threadlocal = (threadlocal_t *)arg;
		if (threadlocal && threadlocal->host) {
				FREE(threadlocal->host);
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
		hostdef_t *host = NULL;

		/* these variables are for the socket */
		struct sockaddr_in addr;
		struct ip_mreq mreq;
		int optval;

		threadlocal_t threadlocal;
		threadlocal.host = NULL;
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

				if ((host = malloc(sizeof(hostdef_t)))==NULL) {
						perror("discover malloc");
						goto cleanup;
				}

				/* this will be deallocated at cleanup */
				threadlocal.host = host;

				/* note that this might be thread cancelation point, but I am not sure */
				if ((nbytes=recvfrom(fd,host,sizeof(hostdef_t),0,(struct sockaddr *)&addr,&addrlen)) < 0) {
						perror("discover recvfrom");
						goto cleanup;
				}

				if (verbose>0) {
						fprintf(stderr, "\n");
						fprintf(stderr, "discover: host->name = %s\n", host->name);
						fprintf(stderr, "discover: host->port = %d\n", host->port);
						fprintf(stderr, "discover: host->id   = %d\n", host->id);
				}

				/* check whether the peer should be listed */
				accept = (accept & ismember_userlist (host->user));
				accept = (accept & ismember_grouplist(host->group));
				accept = (accept & ismember_hostlist (host->name));

				if (!accept) {
						FREE(host);
						continue;
				}

				pthread_mutex_lock(&mutexpeerlist);

				/* remove previous observations of this host from the list */
				/* test the first item on the list */
				if (peerlist)
				{
						found = (peerlist->host->id==host->id);
						if (found) {
								/* delete the first item in the list */
								next = peerlist->next;
								FREE(peerlist);
								peerlist = next;
						}
				}

				/* remove previous observations of this host from the list */
				/* traverse the list */
				peer = peerlist;
				while(peer) {
						/* test the next item on the list */
						next = peer->next;
						if (next) {
								found = (next->host->id==host->id);
								if (found) {
										/* delete the next item in the list */
										peer->next = next->next;
										FREE(next);
										break;
								}
						}
						peer = peer->next;
				}

				/* add the new discovery to the list */
				peer = (peerlist_t *)malloc(sizeof(peerlist_t));
				peer->host = (hostdef_t *)malloc(sizeof(hostdef_t));
				memcpy(peer->host, host, sizeof(hostdef_t));
				FREE(host);

				peer->time      = time(NULL);
				peer->next      = peerlist;
				peerlist        = peer;

				if (verbose>1) {
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
