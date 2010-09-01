#ifndef EXTERN_H
#define EXTERN_H

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

extern pthread_mutex_t mutexhost;
extern hostdef_t *host;

extern pthread_mutex_t mutexpeerlist;
extern peerlist_t *peerlist;

extern pthread_mutex_t mutexjoblist;
extern joblist_t *joblist;

extern pthread_mutex_t mutexuserlist;
extern userlist_t *userlist;

extern pthread_mutex_t mutexgrouplist;
extern grouplist_t *grouplist;

extern pthread_mutex_t mutexhostlist;
extern hostlist_t *hostlist;

extern pthread_mutex_t mutexsmartmem;
extern struct {
		int enabled;
} smartmem;

extern pthread_mutex_t mutexsmartcpu;
extern struct {
		int enabled;
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
		int    n;
		time_t t0;
		int prevhostcount;
		int prevhostid;
		int enabled;
} smartshare;


#endif
