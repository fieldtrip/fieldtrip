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

int smartcpu_info(int *ProcessorCount, float *BogoMips, float *AvgLoad, float *CpuLoad) {
		int verbose = 0;
		void *fp;
		char str[256];
		float AvgLoad1 = 0, AvgLoad5 = 0, AvgLoad15 = 0;

		*ProcessorCount = 0;
		*BogoMips       = 0;
		*AvgLoad        = 0;
		*CpuLoad        = 0;

#if defined (PLATFORM_LINUX)

		/* processor   : 0
		 * vendor_id   : GenuineIntel
		 * cpu family  : 6
		 * model       : 26
		 * model name  : Intel(R) Xeon(R) CPU           W3530  @ 2.80GHz
		 * stepping    : 5
		 * cpu MHz     : 2800.252
		 * cache size  : 8192 KB
		 * physical id : 0
		 * siblings    : 4
		 * core id     : 0
		 * cpu cores   : 4
		 * fpu     : yes
		 * fpu_exception   : yes
		 * cpuid level : 11
		 * wp      : yes
		 * flags       : fpu vme de pse tsc msr pae mce cx8 apic sep mtrr pge mca cmov pat pse36 clflush dts a
		 * bogomips    : 5604.46
		 * clflush size    : 64
		 * cache_alignment : 64
		 * address sizes   : 36 bits physical, 48 bits virtual
		 * power management:
		 */

		if ((fp = fopen("/proc/cpuinfo", "r")) == NULL) {
				fprintf(stderr, "could not open /proc/cpuinfo");
				return -1;
		}

		while (fscanf(fp, "%s :", str) != EOF) {
				if (strncmp(str, "processor", 9)==0) 
						(*ProcessorCount)++;
				if (strncmp(str, "bogomips", 8)==0) 
						fscanf(fp, "%f", BogoMips);
		} /* while */

		fclose(fp);

		/*
		 * 1.85 1.50 1.26 3/221 31765
		 */

		if ((fp = fopen("/proc/loadavg", "r")) == NULL) {
				fprintf(stderr, "could not open /proc/loadavg");
				return -1;
		}

		fscanf(fp, "%f %f %f", &AvgLoad1, &AvgLoad5, &AvgLoad15);
		*AvgLoad = AvgLoad1;

		fclose(fp);

		/*
		 * The columns represent: user nice system idle iowait irq softirq
		 *
		 * cpu  9181729 7423661 2979069 893321421 705481 21805 268125 0
		 * cpu0 2846118 1899136 677048 222954488 90347 0 8216 0
		 * cpu1 2122128 1816939 577506 223890808 63992 1 3941 0
		 * cpu2 2411028 1865444 647992 223263483 254241 8639 24482 0
		 * cpu3 1802453 1842141 1076522 223212640 296899 13163 231485 0
		 * ...
		 */

		if ((fp = fopen("/proc/stat", "r")) == NULL) {
				fprintf(stderr, "could not open /proc/stat");
				return -1;
		}

		int ok, numcpu = 0, elapsed;
		int user, nice, system, idle, iowait, irq, softirq, unknown;
		float load;
		while (!feof(fp)) {
				if (fscanf(fp, "%s ", str) != 1) {
						fprintf(stderr, "smartcpu_info: problem reading from /proc/stat\n");
						continue;
				}

				if (strcmp(str, "cpu")==0) {
						ok = 1;
						ok = ok & (fscanf(fp, "%d", &user    ) == 1);
						ok = ok & (fscanf(fp, "%d", &nice    ) == 1);
						ok = ok & (fscanf(fp, "%d", &system  ) == 1);
						ok = ok & (fscanf(fp, "%d", &idle    ) == 1);
						ok = ok & (fscanf(fp, "%d", &iowait  ) == 1);
						ok = ok & (fscanf(fp, "%d", &irq     ) == 1);
						ok = ok & (fscanf(fp, "%d", &softirq ) == 1);
						ok = ok & (fscanf(fp, "%d", &unknown ) == 1);
						if (!ok) {
								fprintf(stderr, "smartcpu_info: problem reading from /proc/stat\n");
								return -1;
						}
						if (verbose>1)
								fprintf(stderr, "%s %d %d %d %d %d %d %d %d\n", str, user, nice, system, idle, iowait, irq, softirq, unknown);

						pthread_mutex_lock(&mutexprevcpu);
						/* compute the current sum, the previous sum, and the difference */
						elapsed  = user+nice+system+idle+iowait+irq+softirq+unknown;
						elapsed -= prevcpu.user+prevcpu.nice+prevcpu.system+prevcpu.idle+prevcpu.iowait+prevcpu.irq+prevcpu.softirq+prevcpu.unknown;

						if (prevcpu.user==0)
								/* this is the first measurement, the previous values are missing */
								load = 0;
						else
								/* the load is the relative fraction of idle time */
								load = 1.0 -((float)(idle - prevcpu.idle))/((float)elapsed);

						/* update the values that will be used the next time */
						prevcpu.user    = user   ;
						prevcpu.nice    = nice   ;
						prevcpu.system  = system ;
						prevcpu.idle    = idle   ;
						prevcpu.iowait  = iowait ;
						prevcpu.irq     = irq    ;
						prevcpu.softirq = softirq;
						prevcpu.unknown = unknown;
						pthread_mutex_unlock(&mutexprevcpu);

				}
				else if (strncmp(str, "cpu", 3)==0) {
						numcpu++;
						ok = 1;
						ok = ok & (fscanf(fp, "%d", &user    ) == 1);
						ok = ok & (fscanf(fp, "%d", &nice    ) == 1);
						ok = ok & (fscanf(fp, "%d", &system  ) == 1);
						ok = ok & (fscanf(fp, "%d", &idle    ) == 1);
						ok = ok & (fscanf(fp, "%d", &iowait  ) == 1);
						ok = ok & (fscanf(fp, "%d", &irq     ) == 1);
						ok = ok & (fscanf(fp, "%d", &softirq ) == 1);
						ok = ok & (fscanf(fp, "%d", &unknown ) == 1);
						if (!ok) {
								fprintf(stderr, "smartcpu_info: problem reading from /proc/stat\n");
								return -1;
						}
						if (verbose>1)
								fprintf(stderr, "%s %d %d %d %d %d %d %d %d\n", str, user, nice, system, idle, iowait, irq, softirq, unknown);
				}
		} /* while */

		*CpuLoad =  load*numcpu;
		if (verbose>0)
				fprintf(stderr, "smartcpu_info: load = %f %%, numcpu = %d\n", 100*load, numcpu);

		fclose(fp);

		return 0;
#else
		return -1;
#endif

} /* smartcpu_info */

int smartcpu_update(void) {
		peerlist_t *peer;
		int ProcessorCount, NumPeers=0;
		int ok, verbose = 0;
		float BogoMips, AvgLoad, CpuLoad;

		pthread_mutex_lock(&mutexsmartcpu);
		if (!smartcpu.enabled) {
				pthread_mutex_unlock(&mutexsmartcpu);
				return 0;
		}

		pthread_mutex_lock(&mutexhost);
		if (host->status!=STATUS_IDLE & smartcpu.prevstatus!=STATUS_IDLE) {
				pthread_mutex_unlock(&mutexhost);
				pthread_mutex_unlock(&mutexsmartcpu);
				return 0;
		}

		pthread_mutex_lock(&mutexpeerlist);
		/* count the number of idle slave peers on this computer */
		peer = peerlist;
		while(peer) {
				ok = 1;
				ok = ok & (strcmp(peer->ipaddr, "127.0.0.1")==0);
				ok = ok & (peer->host->status == STATUS_IDLE);
				if (ok) {
						NumPeers++;
				}
				peer = peer->next ;
		}
		pthread_mutex_unlock(&mutexpeerlist);

		/* determine the CPU and load details of this computer */
		if (smartcpu_info(&ProcessorCount, &BogoMips, &AvgLoad, &CpuLoad) < 0) {
				pthread_mutex_unlock(&mutexhost);
				return -1;
		}

		/* the numer of idle slaves should not exceed the available free CPUs */
		/* to avoid a race condition with the other idle slaves, the decision is based on two observations */
		if (host->status==STATUS_IDLE & ((float)ProcessorCount-(float)NumPeers-CpuLoad+0.05) < (0-SMARTCPU_TOLERANCE))
				/* increase the evidence to switch from idle to zombie */
				smartcpu.evidence--;
		else if (host->status==STATUS_ZOMBIE & ((float)ProcessorCount-(float)NumPeers-CpuLoad) > (1-SMARTCPU_TOLERANCE)) 
				/* increase the evidence to switch from zombie to idle */
				smartcpu.evidence++;
		else
				/* the current status is fine */
				smartcpu.evidence=0;

		if (smartcpu.evidence<-1) {
				smartcpu.evidence   = 0;
				smartcpu.prevstatus = STATUS_IDLE;
				host->status        = STATUS_ZOMBIE;

				if (verbose>0)
						fprintf(stderr, "smartcpu_update: switching to zombie\n");

				if (verbose>1) {
						fprintf(stderr, "smartcpu_update: ProcessorCount = %d\n", ProcessorCount);
						fprintf(stderr, "smartcpu_update: NumPeers       = %d\n", NumPeers);
						fprintf(stderr, "smartcpu_update: BogoMips       = %.2f\n", BogoMips);
						fprintf(stderr, "smartcpu_update: AvgLoad        = %.2f\n", AvgLoad);
						fprintf(stderr, "smartcpu_update: CpuLoad        = %.2f %%\n", CpuLoad*100);
						fprintf(stderr, "smartcpu_update: host->status   = %d\n", host->status);
				}
		} /* if evidence */

		if (smartcpu.evidence>1) {
				smartcpu.evidence++;
				smartcpu.prevstatus = STATUS_ZOMBIE;
				host->status        = STATUS_IDLE;

				if (verbose>0)
						fprintf(stderr, "smartcpu_update: switching to idle\n");

				if (verbose>1) {
						fprintf(stderr, "smartcpu_update: ProcessorCount = %d\n", ProcessorCount);
						fprintf(stderr, "smartcpu_update: NumPeers       = %d\n", NumPeers);
						fprintf(stderr, "smartcpu_update: BogoMips       = %.2f\n", BogoMips);
						fprintf(stderr, "smartcpu_update: AvgLoad        = %.2f\n", AvgLoad);
						fprintf(stderr, "smartcpu_update: CpuLoad        = %.2f %%\n", CpuLoad*100);
						fprintf(stderr, "smartcpu_update: host->status   = %d\n", host->status);
				}
		} /* if evidence */

		pthread_mutex_unlock(&mutexhost);
		pthread_mutex_unlock(&mutexsmartcpu);

		return 0;
}
