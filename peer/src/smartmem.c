#include <stdio.h>
#include <stdlib.h>

#include "peer.h"
#include "extern.h"
#include "platform_includes.h"

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
				fprintf(stderr, "could not open /proc/meminfo");
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

int smartmem_update(void) {
		peerlist_t* peer;
		unsigned int NumPeers=0;
		UINT64_T MemSuggested=0, MemReserved=0;
		UINT64_T MemTotal=0, MemFree=0, Buffers=0, Cached=0;
		float scale;
		int verbose = 0, status;

		pthread_mutex_lock(&mutexsmartmem);
		if (!smartmem_enabled) {
				pthread_mutex_unlock(&mutexsmartmem);
				return 0;
		}
		pthread_mutex_unlock(&mutexsmartmem);

		/* determine the amount of memory available on this computer */
		if ((status = smartmem_info(&MemTotal, &MemFree, &Buffers, &Cached)) < 0)
				return -1;

		pthread_mutex_lock(&mutexhost);
		pthread_mutex_lock(&mutexpeerlist);

		/* status = 0 means zombie mode, don't accept anything   */
		/* status = 1 means master mode, accept everything       */
		/* status = 2 means idle slave, accept only a single job */
		/* status = 3 means busy slave, don't accept a new job   */
		/* any other status is interpreted as zombie mode        */

		/* determine the amount of memory that is reserved by the other peers on this computer */
		peer = peerlist;
		while(peer) {
				status = 1;
				status = status & (strcmp(peer->ipaddr, "127.0.0.1")==0);
				status = status & (peer->host->id != host->id);
				status = status & (peer->host->status == 2 | peer->host->status == 3); /* include idle and busy slaves */
				if (status) {
						MemReserved += peer->host->memavail;
						NumPeers++;
				}
				peer = peer->next ;
		}

		pthread_mutex_unlock(&mutexhost);
		pthread_mutex_unlock(&mutexpeerlist);

		/* also count this peer */
		NumPeers++;

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

		/* the Buffers and Cached memory are also available to the system */
		MemFree += Buffers;
		MemFree += Cached;

		/* the reserved memory cannot be more than the available memory */
		MemReserved  = (MemReserved > MemFree ? MemFree : MemReserved );

		/* determine the suggested amount of memory for this slave */
		/* this asymptotically approaches the free memory          */
		MemSuggested = ((float)(MemFree - MemReserved)) / (1.0 + scale);

		/* it does not make sense to suggest less than 100MB */
		MemSuggested = (MemSuggested > 1024*1024*100 ? MemSuggested : 1024*1024*100 );

		pthread_mutex_lock(&mutexhost);
		if (host->status==2) {
				host->memavail = MemSuggested;
				if (verbose>0) {
						fprintf(stderr, "NumPeers     = %u\n",   NumPeers);
						fprintf(stderr, "MemFree      = %llu (%f GB)\n", MemFree     , ((float)MemFree     )/(1024*1024*1024));
						fprintf(stderr, "MemReserved  = %llu (%f GB)\n", MemReserved , ((float)MemReserved )/(1024*1024*1024));
						fprintf(stderr, "MemSuggested = %llu (%f GB)\n", MemSuggested, ((float)MemSuggested)/(1024*1024*1024));
				}
		}
		pthread_mutex_unlock(&mutexhost);
		return 0;

}
