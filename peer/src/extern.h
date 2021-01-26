#ifndef EXTERN_H
#define EXTERN_H

extern int syslog_level;

extern pthread_cond_t condstatus;
extern pthread_mutex_t mutexstatus;
extern int udsserverStatus;
extern int tcpserverStatus;
extern int announceStatus;
extern int discoverStatus;
extern int expireStatus;

extern pthread_mutex_t mutexappendcount;
extern int appendcount;

extern pthread_mutex_t mutexsocketcount;
extern int socketcount;

extern pthread_mutex_t mutexthreadcount;
extern int threadcount;

extern pthread_mutex_t mutexconnectioncount;
extern int connectioncount;

extern pthread_mutex_t mutexsmartfreeze;
extern int smartfreeze;

extern pthread_mutex_t mutexhost;
extern hostdef_t *host;

extern pthread_mutex_t mutexpeerlist;
extern peerlist_t *peerlist;

extern pthread_mutex_t mutexjoblist;
extern joblist_t *joblist;

extern pthread_mutex_t mutexallowuserlist;
extern userlist_t *allowuserlist;

extern pthread_mutex_t mutexrefuseuserlist;
extern userlist_t *refuseuserlist;

extern pthread_mutex_t mutexallowgrouplist;
extern grouplist_t *allowgrouplist;

extern pthread_mutex_t mutexrefusegrouplist;
extern grouplist_t *refusegrouplist;

extern pthread_mutex_t mutexallowhostlist;
extern hostlist_t *allowhostlist;

extern pthread_mutex_t mutexrefusehostlist;
extern hostlist_t *refusehostlist;

extern pthread_mutex_t mutexwatchdog;
extern struct {
		int enabled;
		int evidence;
		UINT32_T controllerid;
		time_t time;
		UINT64_T memory;
} watchdog;

extern pthread_mutex_t mutexsmartmem;
extern struct {
		int enabled;
		int freeze;
        UINT64_T memavail;
} smartmem;

extern pthread_mutex_t mutexsmartcpu;
extern struct {
		int enabled;
		int freeze;
		int prevstatus;
		int evidence;
} smartcpu;

extern pthread_mutex_t mutexprevcpu;
extern struct {
		int user, nice, system, idle, iowait, irq, softirq, unknown;
} prevcpu;

extern pthread_mutex_t mutexsmartshare;
extern smartsharelist_t *smartsharelist;
extern struct {
		int enabled;
		int prevhostcount;
		int prevhostid;
		int n;
		time_t time;
} smartshare;


#endif
