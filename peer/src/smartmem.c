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

#include "peer.h"
#include "extern.h"
#include "platform_includes.h"

/* return value 0 if ok, -1 if error or unsupported platform */
int smartmem_info(UINT64_T *MemTotal, UINT64_T *MemFree, UINT64_T *Buffers, UINT64_T *Cached) {
		void *fp;
		char str[256];

		*MemTotal = UINT32_MAX;
		*MemFree  = 0;
		*Buffers  = 0;
		*Cached   = 0;

#if defined (PLATFORM_LINUX)

		/* get the memory information from /proc/meminfo
		 *
		 * MemTotal:     32949684 kB
		 * MemFree:       4200556 kB
		 * Buffers:         73612 kB
		 * Cached:       26112456 kB
		 * SwapCached:     104324 kB
		 * Active:       10769324 kB
		 * Inactive:     17096716 kB
		 * HighTotal:           0 kB
		 * HighFree:            0 kB
		 * LowTotal:     32949684 kB
		 * LowFree:       4200556 kB
		 * SwapTotal:    62857388 kB
		 * SwapFree:     62701340 kB
		 * Dirty:         3110792 kB
		 * Writeback:           0 kB
		 * AnonPages:     1677744 kB
		 * Mapped:         101508 kB
		 * Slab:           788800 kB
		 * PageTables:      26740 kB
		 * NFS_Unstable:        0 kB
		 * Bounce:              0 kB
		 * CommitLimit:  79332228 kB
		 * Committed_AS:  2660040 kB
		 * VmallocTotal: 34359738367 kB
		 * VmallocUsed:    293732 kB
		 * VmallocChunk: 34359443995 kB
		 * HugePages_Total:     0
		 * HugePages_Free:      0
		 * HugePages_Rsvd:      0
		 * Hugepagesize:     2048 kB
		 */

		if ((fp = fopen("/proc/meminfo", "r")) == NULL) {
				DEBUG(LOG_ERR, "smartmem_info: could not open /proc/meminfo");
				return -1;
		}

		while (fscanf(fp, "%s", str) != EOF) {
				if (strcmp(str, "MemTotal:")==0) 
						fscanf(fp, "%llu", MemTotal);
				if (strcmp(str, "MemFree:")==0) 
						fscanf(fp, "%llu", MemFree);
				if (strcmp(str, "Buffers:")==0) 
						fscanf(fp, "%llu", Buffers);
				if (strcmp(str, "Cached:")==0) 
						fscanf(fp, "%llu", Cached);

		} /* while */

		fclose(fp);

		/* convert from kB into bytes */
		(*MemTotal) *= 1024;
		(*MemFree)  *= 1024;
		(*Buffers)  *= 1024;
		(*Cached)   *= 1024;

		return 0;
#else
		return -1;
#endif
} /* smartmem_info */

/* return value 0 if ok, -1 if error or unsupported platform */
int smartmem_update(void) {
		peerlist_t* peer;
		unsigned int NumPeers=0;
		UINT64_T MemSuggested=0, MemReserved=0;
		UINT64_T MemTotal=0, MemFree=0, Buffers=0, Cached=0;
		float scale;
		int ok;

		pthread_mutex_lock(&mutexsmartmem);
		if (!smartmem.enabled || smartmem.freeze) {
				/* don't update if smartmem is disabled */
				/* the freeze flag is enabled during the processing of an incoming job in tcpsocket */
				pthread_mutex_unlock(&mutexsmartmem);
				return 0;
		}

		pthread_mutex_lock(&mutexhost);
		if (host->status!=STATUS_IDLE) {
				/* don't update if the status is something else than IDLE */
				pthread_mutex_unlock(&mutexhost);
				pthread_mutex_unlock(&mutexsmartmem);
				return 0;
		}

		/* determine the amount of memory available on this computer */
		if (smartmem_info(&MemTotal, &MemFree, &Buffers, &Cached) < 0) {
				pthread_mutex_unlock(&mutexhost);
				pthread_mutex_unlock(&mutexsmartmem);
				return -1;
		}

		DEBUG(LOG_INFO, "smartmem: host->memavail = %llu", host->memavail);
		DEBUG(LOG_DEBUG, "smartmem: MemTotal       = %llu (%f GB)", MemTotal    , ((float)MemTotal    )/(1024*1024*1024));
		DEBUG(LOG_DEBUG, "smartmem: MemFree        = %llu (%f GB)", MemFree     , ((float)MemFree     )/(1024*1024*1024));
		DEBUG(LOG_DEBUG, "smartmem: Buffers        = %llu (%f GB)", Buffers     , ((float)Buffers     )/(1024*1024*1024));
		DEBUG(LOG_DEBUG, "smartmem: Cached         = %llu (%f GB)", Cached      , ((float)Cached      )/(1024*1024*1024));

		pthread_mutex_lock(&mutexpeerlist);
		/* determine the amount of memory that is reserved by the idle workers on this computer */
		peer = peerlist;
		while(peer) {
				ok = 1;
				ok = ok && (strcmp(peer->ipaddr, "127.0.0.1")==0);
				ok = ok && (peer->host->id != host->id);
				ok = ok && (peer->host->status == STATUS_IDLE);
				if (ok) {
						MemReserved += peer->host->memavail;
						NumPeers++;
				}
				peer = peer->next ;
		}
		pthread_mutex_unlock(&mutexpeerlist);

		/* include this peer in the count */
		NumPeers++;

		/* the Buffers and Cached memory are also available to the system */
		MemFree += Buffers;
		MemFree += Cached;

		/* determine the scale of the memory update (for iterative improvements) */
		scale = ((float)host->memavail) / ((float)MemFree);

		/* apply a heuristic adjustment to improve total memory usage (in case of few peers) and to speed up convergence (in case of many peers) */
		switch (NumPeers)
		{
				case 1:
						scale = 0.01 * scale;
						break;
				case 2:
						scale = 0.1 * scale;
						break;
				case 3:
						scale = 0.3 * scale;
						break;
				case 4:
						scale = 0.5 * scale;
						break;
				default:
						scale = 1.0 * scale;
						break;
		} /* switch */

		/* the reserved memory by the other idle workers should not be more than the available memory */
		MemReserved  = (MemReserved > MemFree ? MemFree : MemReserved );

		/* determine the suggested amount of memory for this worker */
		/* this asymptotically approaches the free memory          */
		MemSuggested = ((float)(MemFree - MemReserved)) / (1.0 + scale);

		/* it does not make sense to suggest less than a certain minimum amount */
		MemSuggested = (MemSuggested > SMARTMEM_MINIMUM ? MemSuggested : SMARTMEM_MINIMUM );

		/* don't exceed the user-specified maximum allowed amount */
		MemSuggested = (MemSuggested < smartmem.memavail ? MemSuggested : smartmem.memavail );

		DEBUG(LOG_DEBUG, "smartmem: NumPeers       = %u",   NumPeers);
		DEBUG(LOG_DEBUG, "smartmem: MemReserved    = %llu (%f GB)", MemReserved , ((float)MemReserved )/(1024*1024*1024));
		DEBUG(LOG_DEBUG, "smartmem: MemSuggested   = %llu (%f GB)", MemSuggested, ((float)MemSuggested)/(1024*1024*1024));

		host->memavail = MemSuggested;

		pthread_mutex_unlock(&mutexhost);
		pthread_mutex_unlock(&mutexsmartmem);

		return 0;
} /* smartmem_update */
