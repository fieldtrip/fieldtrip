/*
 * Copyright (C) 2008-2010, Robert Oostenveld
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

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>       /* for strerror */

#include "peer.h"
#include "extern.h"
#include "platform_includes.h"

int bufread(int s, void *buf, int numel) {
		int numcall = 0, numthis = 0, numread = 0;

		while (numread<numel) {

				numthis = recv(s, (char*)buf+numread, numel-numread, 0);
				if (numthis<0) {
						perror("bufread");
						DEBUG(LOG_ERR, "error: bufread");
						break;
				}
				else if (numthis == 0)
						break;

				DEBUG(LOG_DEBUG, "bufread: read %d bytes", numthis);
				numread += numthis;
				numcall++;
				usleep(1000);
		}
		DEBUG(LOG_DEBUG, "bufread: reading the complete buffer required %d calls", numcall);
		return numread;
}

int bufwrite(int s, void *buf, int numel) {
		int numcall = 0, numthis = 0, numwrite = 0;

		DEBUG(LOG_DEBUG, "bufwrite: request for %d bytes", numel);

		while (numwrite<numel) {

				numthis = send(s, (char*)buf+numwrite, numel-numwrite, 0);
				if (numthis<0) {
						perror("bufwrite");
						DEBUG(LOG_ERR, "error: bufwrite");
						break;
				}
				else if(numthis == 0)
						break;

				DEBUG(LOG_DEBUG, "bufwrite: wrote %d bytes", numthis);
				numwrite += numthis;
				numcall++;
				usleep(1000);
		}
		DEBUG(LOG_DEBUG, "bufwrite: writing the complete buffer required %d calls", numcall);
		return numwrite;
}

int append(void **buf1, int bufsize1, void *buf2, int bufsize2) {

		pthread_mutex_lock(&mutexappendcount);
		appendcount++;
		DEBUG(LOG_DEBUG, "append: appendcount = %d", appendcount);
		pthread_mutex_unlock(&mutexappendcount);

		if (((*buf1)!=NULL) && (bufsize1==0)) {
				perror("append err1");
				DEBUG(LOG_ERR, "error: append err1");
				return -1;
		}
		else if (((*buf1)==NULL) && (bufsize1!=0)) {
				perror("append err2");
				DEBUG(LOG_ERR, "error: append err2");
				return -1;
		}

		if ((*buf1)==NULL) {
				DEBUG(LOG_DEBUG, "append: allocating %d bytes", bufsize2);
				(*buf1) = malloc(bufsize2);
		}
		else if ((*buf1)!=NULL) {
				DEBUG(LOG_DEBUG, "append: reallocating from %d to %d bytes", bufsize1, bufsize1+bufsize2);
				(*buf1) = realloc((*buf1), bufsize1+bufsize2);
		}

		memcpy((char*)(*buf1)+bufsize1, buf2, bufsize2);
		return (bufsize1+bufsize2);
}

int jobcount(void) {
		int jobcount = 0;
		joblist_t *job;
		pthread_mutex_lock(&mutexjoblist);
		job = joblist;
		while (job) {
				jobcount++;
				job = job->next;
		}
		pthread_mutex_unlock(&mutexjoblist);
		return jobcount;
}

int peercount(void) {
		int peercount = 0;
		peerlist_t *peer;
		pthread_mutex_lock(&mutexpeerlist);
		peer = peerlist;
		while (peer) {
				peercount++;
				peer = peer->next;
		}
		pthread_mutex_unlock(&mutexpeerlist);
		return peercount;
}

int hoststatus(void) {
		int status = -1;
		pthread_mutex_lock(&mutexhost);
		if (host)
				status = host->status;
		pthread_mutex_unlock(&mutexhost);
		return status;
}

void clear_peerlist(void) {
		peerlist_t *peer = NULL;
		pthread_mutex_lock(&mutexpeerlist);
		peer = peerlist;
		while (peer) {
				peerlist = peer->next;
				FREE(peer->host);
				FREE(peer);
				peer = peerlist;
		}
		pthread_mutex_unlock(&mutexpeerlist);
}

void clear_joblist(void) {
		joblist_t *job = NULL;
		pthread_mutex_lock(&mutexjoblist);
		job = joblist;
		while (job) {
				joblist = job->next;
				FREE(job->job);
				FREE(job->host);
				FREE(job->arg);
				FREE(job->opt);
				FREE(job);
				job = joblist;
		}
		pthread_mutex_unlock(&mutexjoblist);
}

void clear_userlist(void) {
		userlist_t *user = NULL;
		pthread_mutex_lock(&mutexuserlist);
		user = userlist;
		while (user) {
				userlist = user->next;
				FREE(user->name);
				FREE(user);
				user = userlist;
		}
		pthread_mutex_unlock(&mutexuserlist);
}

void clear_grouplist(void) {
		grouplist_t *group = NULL;
		pthread_mutex_lock(&mutexgrouplist);
		group = grouplist;
		while (group) {
				grouplist = group->next;
				FREE(group->name);
				FREE(group);
				group = grouplist;
		}
		pthread_mutex_unlock(&mutexgrouplist);
}

void clear_hostlist(void) {
		hostlist_t *listitem = NULL;
		pthread_mutex_lock(&mutexhostlist);
		listitem = hostlist;
		while (listitem) {
				hostlist = listitem->next;
				FREE(listitem->name);
				FREE(listitem);
				listitem = hostlist;
		}
		pthread_mutex_unlock(&mutexhostlist);
}

void clear_smartsharelist(void) {
		smartsharelist_t *listitem = NULL;
		pthread_mutex_lock(&mutexsmartshare);
		listitem = smartsharelist;
		while (listitem) {
				smartsharelist = listitem->next;
				FREE(listitem);
				listitem = smartsharelist;
		}
		pthread_mutex_unlock(&mutexsmartshare);
}

void check_datatypes() {
		/* check datatypes */
		if (WORDSIZE_CHAR    !=1) { PANIC("invalid size of CHAR    (%d)", WORDSIZE_CHAR   );  }
		if (WORDSIZE_UINT8   !=1) { PANIC("invalid size of UINT8   (%d)", WORDSIZE_UINT8  );  }
		if (WORDSIZE_UINT16  !=2) { PANIC("invalid size of UINT16  (%d)", WORDSIZE_UINT16 );  }
		if (WORDSIZE_UINT32  !=4) { PANIC("invalid size of UINT32  (%d)", WORDSIZE_UINT32 );  }
		if (WORDSIZE_UINT64  !=8) { PANIC("invalid size of UINT64  (%d)", WORDSIZE_UINT64 );  }
		if (WORDSIZE_INT8    !=1) { PANIC("invalid size of INT8    (%d)", WORDSIZE_INT8   );  }
		if (WORDSIZE_INT16   !=2) { PANIC("invalid size of INT16   (%d)", WORDSIZE_INT16  );  }
		if (WORDSIZE_INT32   !=4) { PANIC("invalid size of INT32   (%d)", WORDSIZE_INT32  );  }
		if (WORDSIZE_INT64   !=8) { PANIC("invalid size of INT64   (%d)", WORDSIZE_INT64  );  }
		if (WORDSIZE_FLOAT32 !=4) { PANIC("invalid size of FLOAT32 (%d)", WORDSIZE_FLOAT32);  }
		if (WORDSIZE_FLOAT64 !=8) { PANIC("invalid size of FLOAT64 (%d)", WORDSIZE_FLOAT64);  }
		if (sizeof(hostdef_t)!=976) { PANIC("invalid size of hostdef_t (%d)", sizeof(hostdef_t) );  }
}

