#ifndef EXTERN_H
#define EXTERN_H

extern pthread_mutex_t mutexstatus;
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

extern pthread_mutex_t mutexfairshare;
extern struct {
  float a;
  float b;
  float c;
  float d;
  int   n;
  int prevhostcount;
  int prevhostid;
  time_t t0;
} param;

#endif
