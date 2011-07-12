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

/* the following is for getmem */
#if defined (PLATFORM_OSX)
#include <mach/task.h>
#elif defined(PLATFORM_WIN32) || defined(PLATFORM_WIN64)
#include <windows.h>
#include <psapi.h>
#elif defined(PLATFORM_LINUX)
#endif

int threadsleep(float t) {
	#ifdef WIN32
		Sleep((int) (t*1000.0));
		return 0;
	#else
		int retval = 0;
		struct timespec req, rem;
		/* split in seconds and nanoseconds */
		req.tv_sec  = (int)t;
		req.tv_nsec = (int)(1000000000.0 * (t - (int)t));
		retval = nanosleep(&req, &rem);
		return retval;
	#endif
}

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
				threadsleep(0.001);
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
				threadsleep(0.001);
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

void clear_allowuserlist(void) {
		userlist_t *user = NULL;
		pthread_mutex_lock(&mutexallowuserlist);
		user = allowuserlist;
		while (user) {
				allowuserlist = user->next;
				FREE(user->name);
				FREE(user);
				user = allowuserlist;
		}
		pthread_mutex_unlock(&mutexallowuserlist);
}

void clear_allowgrouplist(void) {
		grouplist_t *group = NULL;
		pthread_mutex_lock(&mutexallowgrouplist);
		group = allowgrouplist;
		while (group) {
				allowgrouplist = group->next;
				FREE(group->name);
				FREE(group);
				group = allowgrouplist;
		}
		pthread_mutex_unlock(&mutexallowgrouplist);
}

void clear_allowhostlist(void) {
		hostlist_t *listitem = NULL;
		pthread_mutex_lock(&mutexallowhostlist);
		listitem = allowhostlist;
		while (listitem) {
				allowhostlist = listitem->next;
				FREE(listitem->name);
				FREE(listitem);
				listitem = allowhostlist;
		}
		pthread_mutex_unlock(&mutexallowhostlist);
}

void clear_refuseuserlist(void) {
		userlist_t *user = NULL;
		pthread_mutex_lock(&mutexrefuseuserlist);
		user = refuseuserlist;
		while (user) {
				refuseuserlist = user->next;
				FREE(user->name);
				FREE(user);
				user = refuseuserlist;
		}
		pthread_mutex_unlock(&mutexrefuseuserlist);
}

void clear_refusegrouplist(void) {
		grouplist_t *group = NULL;
		pthread_mutex_lock(&mutexrefusegrouplist);
		group = refusegrouplist;
		while (group) {
				refusegrouplist = group->next;
				FREE(group->name);
				FREE(group);
				group = refusegrouplist;
		}
		pthread_mutex_unlock(&mutexrefusegrouplist);
}

void clear_refusehostlist(void) {
		hostlist_t *listitem = NULL;
		pthread_mutex_lock(&mutexrefusehostlist);
		listitem = refusehostlist;
		while (listitem) {
				refusehostlist = listitem->next;
				FREE(listitem->name);
				FREE(listitem);
				listitem = refusehostlist;
		}
		pthread_mutex_unlock(&mutexrefusehostlist);
}

void check_datatypes() {
		/* check datatypes */
		if (WORDSIZE_CHAR    !=1) { PANIC("invalid size of CHAR    (%lu)", WORDSIZE_CHAR   );  }
		if (WORDSIZE_UINT8   !=1) { PANIC("invalid size of UINT8   (%lu)", WORDSIZE_UINT8  );  }
		if (WORDSIZE_UINT16  !=2) { PANIC("invalid size of UINT16  (%lu)", WORDSIZE_UINT16 );  }
		if (WORDSIZE_UINT32  !=4) { PANIC("invalid size of UINT32  (%lu)", WORDSIZE_UINT32 );  }
		if (WORDSIZE_UINT64  !=8) { PANIC("invalid size of UINT64  (%lu)", WORDSIZE_UINT64 );  }
		if (WORDSIZE_INT8    !=1) { PANIC("invalid size of INT8    (%lu)", WORDSIZE_INT8   );  }
		if (WORDSIZE_INT16   !=2) { PANIC("invalid size of INT16   (%lu)", WORDSIZE_INT16  );  }
		if (WORDSIZE_INT32   !=4) { PANIC("invalid size of INT32   (%lu)", WORDSIZE_INT32  );  }
		if (WORDSIZE_INT64   !=8) { PANIC("invalid size of INT64   (%lu)", WORDSIZE_INT64  );  }
		if (WORDSIZE_FLOAT32 !=4) { PANIC("invalid size of FLOAT32 (%lu)", WORDSIZE_FLOAT32);  }
		if (WORDSIZE_FLOAT64 !=8) { PANIC("invalid size of FLOAT64 (%lu)", WORDSIZE_FLOAT64);  }
		if (sizeof(current_t)!=(416)) { PANIC("invalid size of current_t (%lu)", sizeof(current_t) );  }
		if (sizeof(hostdef_t)!=(552+416)) { PANIC("invalid size of hostdef_t (%lu)", sizeof(hostdef_t) );  }
}

#if defined (PLATFORM_OSX)
int getmem (uint64_t *rss, uint64_t *vs) {
		task_t task = MACH_PORT_NULL;
		struct task_basic_info t_info;
		mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;
		if (KERN_SUCCESS != task_info(mach_task_self(), TASK_BASIC_INFO, (task_info_t)&t_info, &t_info_count)) {
				return -1;
		}
		*rss = t_info.resident_size;
		*vs  = t_info.virtual_size;
		return 0;
}
#elif defined(PLATFORM_WIN32) || defined(PLATFORM_WIN64)
int getmem (uint64_t *rss, uint64_t *vs) {
		/* no idea how to get the memory information on a windows computer */
		*rss = 0;
		*vs  = 0;
		return -1;
}
#elif defined(PLATFORM_LINUX)
int getmem (uint64_t *rss, uint64_t *vs) {
		FILE *fp;
		uint64_t val1 = 0, val2 = 0;
		if ((fp = fopen("/proc/self/statm", "r")) == NULL) {
				DEBUG(LOG_ERR, "could not open /proc/self/statm");
				return -1;
		}
		/* read the information from /proc/self/statm
		   size       total program size
		   resident   resident set size
		   share      shared pages
		   text       text (code)
		   lib        library
		   data       data/stack
		   dt         dirty pages (unused in Linux 2.6)
		 */
		if (fscanf(fp, "%llu%llu", &val1, &val2 )!=2) {
				DEBUG(LOG_WARNING, "could not read all elements from /proc/self/statm");
				val1 = 0;
				val2 = 0;
		}
		/* these seem to be in 4096 byte blocks */
		*vs  = val1 * 4096; 
		*rss = val2 * 4096; 
		fclose(fp);
		return 0;
}
#endif

